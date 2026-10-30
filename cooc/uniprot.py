#!/usr/bin/env python3
"""
Module: uniprot
Purpose: Fetch a UniProt reference proteome (one canonical protein per gene) with its Pfam domains.
Author: Dr. Arli Aditya Parikesit
Date: 2026
Parameters:
    - upid: UniProt reference proteome ID, e.g. UP000005640 (Homo sapiens)
Notes:
    - Genes come from the proteome's *.gene2acc file; the first accession of each gene group is the
      canonical protein. Proteins with an unknown gene ('-') are excluded from both n and the annotation.
    - Domains are UniProtKB Pfam cross-references, named via Pfam-A.clans.tsv (Pfam ID + description),
      i.e. database annotation ("SF" in the 2011 paper), not de novo prediction.
References:
    - The UniProt Consortium (2025). Nucleic Acids Res 53:D609-D617. doi:10.1093/nar/gkae1010
    - Blum M et al. (2025). InterPro: the protein sequence classification resource in 2025.
      Nucleic Acids Res 53:D444-D456. doi:10.1093/nar/gkae1082
"""
import gzip
import json
import re
import time
import urllib.parse

import pandas as pd
import requests

UPID_RE = re.compile(r"UP\d{9}")
REST = "https://rest.uniprot.org"
FTP = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes"
PFAM_CLANS = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz"
RETRY_WAIT = 10  # seconds, multiplied by the attempt number


def _get(url, timeout=900, tries=3):
    """GET and gunzip; retry when the connection drops or a gzip body arrives truncated (EOFError).
    Large UniProt streams run for minutes and are sometimes cut off by servers or proxies."""
    for attempt in range(1, tries + 1):
        try:
            r = requests.get(url, timeout=timeout)
            r.raise_for_status()
            data = r.content
            return (gzip.decompress(data) if data[:2] == b"\x1f\x8b" else data).decode()
        except (requests.ConnectionError, requests.Timeout, requests.exceptions.ChunkedEncodingError,
                EOFError, gzip.BadGzipFile) as e:
            if attempt == tries:
                raise RuntimeError(f"Download failed after {tries} attempts ({e}): {url.split('?')[0]}") from e
            time.sleep(RETRY_WAIT * attempt)


def canonical_genes(gene2acc):
    """gene2acc text -> {canonical accession: gene group}. The group is column 3 (the canonical
    accession's gene); column 1 can differ, e.g. for read-through loci. First row of a group is canonical."""
    out, seen = {}, set()
    for line in gene2acc.splitlines():
        f = line.split("\t")
        if len(f) >= 3 and f[2] != "-" and f[2] not in seen:
            seen.add(f[2])
            out[f[1]] = f[2]
    return out


def pfam_names(clans_tsv):
    """Pfam-A.clans.tsv (acc, clan acc, clan id, Pfam id, description) -> {acc: 'id acc description'}."""
    names = {}
    for line in clans_tsv.splitlines():
        f = line.split("\t")
        if len(f) >= 5:
            names[f[0]] = f"{f[3]} {f[0]} {f[4]}"
    return names


def pfam_annotation(stream_tsv, canon, names):
    """UniProt TSV (Entry, Pfam) -> DataFrame[gene_id, domain] for canonical proteins only."""
    rows = []
    for line in stream_tsv.splitlines()[1:]:
        acc, _, pfams = line.partition("\t")
        if acc in canon:
            rows.extend((canon[acc], names.get(pf, pf)) for pf in pfams.split(";") if pf)
    return pd.DataFrame(rows, columns=["gene_id", "domain"]).drop_duplicates()


def fetch_reference_proteome(upid):
    """Download one reference proteome; return a dataset dict for cooc.core.analyze."""
    upid = upid.strip().upper()
    if not UPID_RE.fullmatch(upid):
        raise ValueError("Proteome ID must look like UP000005640")
    meta = json.loads(_get(f"{REST}/proteomes/{upid}?format=json"))
    taxid, division = meta["taxonomy"]["taxonId"], meta["superkingdom"].capitalize()
    canon = canonical_genes(_get(f"{FTP}/{division}/{upid}/{upid}_{taxid}.gene2acc.gz"))
    if not canon:
        raise ValueError(f"No gene mapping for {upid}")
    names = pfam_names(_get(PFAM_CLANS))
    query = urllib.parse.urlencode(dict(query=f"proteome:{upid} AND database:pfam",
                                        fields="accession,xref_pfam", format="tsv", compressed="true"))
    annot = pfam_annotation(_get(f"{REST}/uniprotkb/stream?{query}"), canon, names)
    return dict(species=f"{meta['taxonomy']['scientificName']} ({upid})", source="UniProt-Pfam",
                n=len(canon), annot=annot)
