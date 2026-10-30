#!/usr/bin/env python3
"""
Module: pipeline
Purpose: De novo gene prediction (AUGUSTUS) and HMM domain annotation (pyhmmer) of predicted proteins.
Author: Dr. Arli Aditya Parikesit
Date: 2026
Parameters:
    - species: AUGUSTUS parameter set (e.g. human, fly, arabidopsis); replaces GENSCAN of the 2011 paper
    - evalue: per-sequence inclusion E-value for hmmsearch
    - gathering: use Pfam gathering (GA) bit-score cutoffs instead of E-value
References:
    - Stanke M, Morgenstern B (2005). Nucleic Acids Res 33:W465-W467. doi:10.1093/nar/gki458
    - Larralde M, Zeller G (2023). Bioinformatics 39(5):btad214. doi:10.1093/bioinformatics/btad214
"""
import os
import re
import shutil
import subprocess
import sys
import threading

import pandas as pd

SPECIES_RE = re.compile(r"[A-Za-z0-9_]{1,64}")
# Binary comes from the environment only, never from the web UI.
AUGUSTUS_BIN = os.environ.get("AUGUSTUS_BIN", "augustus")
# ponytail: one pipeline run per server process; use a job queue if several users need parallel runs.
GP_LOCK = threading.Lock()


def augustus_path():
    return shutil.which(AUGUSTUS_BIN)


def run_augustus(genome_path, species, timeout=1800):
    """Run AUGUSTUS on a genome FASTA and return its GFF3 text (with protein sequences)."""
    if not SPECIES_RE.fullmatch(species):
        raise ValueError("AUGUSTUS species must be alphanumeric/underscore")
    exe = augustus_path()
    if exe is None:
        raise FileNotFoundError(f"AUGUSTUS binary '{AUGUSTUS_BIN}' not found; set AUGUSTUS_BIN")
    proc = subprocess.run([exe, f"--species={species}", "--protein=on", "--gff3=on", str(genome_path)],
                          capture_output=True, text=True, timeout=timeout)
    if proc.returncode != 0:
        print(proc.stderr[-2000:], file=sys.stderr)  # details stay in the server log
        raise RuntimeError(f"AUGUSTUS failed with exit code {proc.returncode}; see server log")
    return proc.stdout


def augustus_proteins(gff_text):
    """Extract '# protein sequence = [...]' blocks, keyed by the preceding transcript ID."""
    prots, tid, buf = {}, None, None
    for line in gff_text.splitlines():
        if buf is not None:
            buf.append(line.lstrip("# ").rstrip("]"))
            if line.rstrip().endswith("]"):
                prots[tid], buf = "".join(buf), None
            continue
        f = line.split("\t")
        if len(f) == 9 and f[2] in ("mRNA", "transcript"):
            m = re.search(r"ID=([^;]+)", f[8])
            tid = m.group(1) if m else f[8].strip()
        elif line.startswith("# protein sequence = ["):
            seq = line.split("[", 1)[1].rstrip()
            if tid is None:
                tid = f"t{len(prots) + 1}"
            if seq.endswith("]"):
                prots[tid] = seq[:-1]
            else:
                buf = [seq]
    return prots


def hmm_annotate(proteins, hmm_path, evalue=1e-3, gathering=False, cpus=0):
    """hmmsearch every profile in hmm_path against proteins; return DataFrame[gene_id, domain]."""
    from pyhmmer import easel, hmmer, plan7

    alphabet = easel.Alphabet.amino()
    seqs = [easel.TextSequence(name=k, sequence=v).digitize(alphabet) for k, v in proteins.items()]
    opts = {"bit_cutoffs": "gathering"} if gathering else {"E": evalue, "incE": evalue}
    rows = []
    with plan7.HMMFile(hmm_path) as hmms:
        for hits in hmmer.hmmsearch(hmms, seqs, cpus=cpus, **opts):
            q = hits.query
            label = " ".join(x for x in (q.name, q.accession, q.description) if x)
            rows.extend((hit.name, label) for hit in hits.included)
    return pd.DataFrame(rows, columns=["gene_id", "domain"])
