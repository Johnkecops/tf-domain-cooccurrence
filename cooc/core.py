#!/usr/bin/env python3
"""
Module: core
Purpose: Parse domain annotations, assign TF domain classes, and test pairwise co-occurrence.
Author: Dr. Arli Aditya Parikesit
Date: 2026
Parameters:
    - n: total number of genes (proteins) in the genome, including those without domains
    - alpha: significance level applied to Benjamini-Hochberg q-values
    - min_e: expectation below which avoidance is reported as undetectable
References:
    - Parikesit AA, Prohaska SJ, Stadler PF (2011). Quantitative Measurement of Genome-wide
      Protein Domain Co-occurrence of Transcription Factors. Leipzig Research Festival for Life Science.
"""
import io
from itertools import combinations

import numpy as np
import pandas as pd
import regex as re
from scipy.stats import false_discovery_control, poisson

# Ordered: the first matching class wins, so specific zinc-binding folds (KRAB, RING, PHD)
# are not also counted as generic znf and do not create artificial znf co-occurrences.
DEFAULT_CLASSES = [
    ("krab", r"KRAB"),
    ("ring", r"RING|U-box|C3HC4"),
    ("phd", r"PHD"),
    ("bromo", r"(?i)bromo"),
    ("brct", r"BRCT"),
    ("hmg", r"HMG[_-]box|(?i:high mobility group)"),
    ("bzip", r"bZIP|(?i:leucine[ -]zipper)"),
    ("wing", r"(?i)winged[ -]helix|wHTH"),
    ("wd40", r"WD40|(?i:WD[ -]?repeat)"),
    ("tpr", r"TPR|(?i:tetratricopeptide)"),
    ("znf", r"zf-|(?i:zinc[ -]finger)|C2H2"),
]

MAX_REGEX_LEN = 200
# User-editable patterns: a per-match timeout stops catastrophic backtracking from stalling the server.
REGEX_TIMEOUT = 0.05


def compile_classes(classes):
    """[(name, regex)] -> [(name, compiled)]; raises ValueError on bad input."""
    out = []
    for name, pattern in classes:
        name, pattern = str(name).strip(), str(pattern).strip()
        if not name or not pattern:
            continue
        if len(pattern) > MAX_REGEX_LEN:
            raise ValueError(f"Regex for '{name}' longer than {MAX_REGEX_LEN} characters")
        try:
            out.append((name, re.compile(pattern)))
        except re.error as e:
            raise ValueError(f"Invalid regex for '{name}': {e}") from e
    return out


def _decode(data):
    return data.decode("utf-8", errors="replace") if isinstance(data, bytes) else data


def read_annotation(data, filename=""):
    """Parse InterProScan TSV, HMMER hmmscan --domtblout, or a 2-column gene/domain table.

    Returns DataFrame[gene_id, domain]."""
    text = _decode(data)
    lines = [l for l in text.splitlines() if l.strip() and not l.startswith("#")]
    if not lines:
        raise ValueError("Annotation file is empty")
    first = lines[0]
    if filename.endswith((".domtbl", ".domtblout")) or len(first.split()) >= 23 and "\t" not in first:
        # hmmscan: target (domain) = col 0/1, query (protein) = col 3, description = col 22+
        rows = []
        for l in lines:
            f = l.split()
            if len(f) < 23:
                raise ValueError("Malformed domtblout line: " + l[:80])
            rows.append((f[3], " ".join([f[0], f[1], " ".join(f[22:])])))
        return pd.DataFrame(rows, columns=["gene_id", "domain"])
    if first.count("\t") >= 10:
        # InterProScan TSV: protein, md5, len, analysis, sig acc, sig desc, ..., col 12 InterPro desc
        rows = []
        for l in lines:
            f = l.split("\t")
            rows.append((f[0], " ".join(x for x in f[4:6] + f[12:13] if x and x != "-")))
        return pd.DataFrame(rows, columns=["gene_id", "domain"])
    header = 0 if "gene" in first.lower() else None
    df = pd.read_csv(io.StringIO(text), sep=None, engine="python", header=header, comment="#", dtype=str)
    if df.shape[1] < 2:
        raise ValueError("Generic table needs two columns: gene_id, domain")
    df = df.iloc[:, :2]
    df.columns = ["gene_id", "domain"]
    return df.dropna()


def count_fasta(data):
    return sum(1 for l in _decode(data).splitlines() if l.startswith(">"))


def classify(annot, classes):
    """Assign each domain hit to its first matching class; return unique (gene_id, cls)."""
    compiled = compile_classes(classes)

    def first(domain):
        try:
            return next((name for name, rx in compiled if rx.search(domain, timeout=REGEX_TIMEOUT)), None)
        except TimeoutError as e:
            raise ValueError(f"Class regex timed out on domain '{domain[:60]}'; simplify the pattern") from e

    cls = annot["domain"].astype(str).map(first)
    out = pd.DataFrame({"gene_id": annot["gene_id"].astype(str), "cls": cls}).dropna()
    return out.drop_duplicates().reset_index(drop=True)


def pair_stats(gene_cls, n):
    """X, Y, F, E = XY/n and one-sided Poisson P for every pair of classes present."""
    genes = gene_cls["gene_id"].nunique()
    if n < genes:
        raise ValueError(f"n ({n}) is smaller than the number of annotated genes ({genes})")
    sets = gene_cls.groupby("cls")["gene_id"].agg(set)
    rows = []
    for a, b in combinations(sorted(sets.index), 2):
        X, Y = len(sets[a]), len(sets[b])
        F = len(sets[a] & sets[b])
        E = X * Y / n
        if F < E:
            direction, p = "avoidance", poisson.cdf(F, E)
        else:
            direction, p = "co-occurrence", poisson.sf(F - 1, E)
        rows.append(dict(pair=f"{a}-{b}", x=a, y=b, X=X, Y=Y, n=n, E=E, F=F,
                         log2_F_over_E=np.log2((F + 0.5) / (E + 0.5)), direction=direction, p=p))
    return pd.DataFrame(rows)


def analyze(datasets, classes=DEFAULT_CLASSES, alpha=0.05, min_e=1.0):
    """datasets: iterable of dicts with species, source, n, annot. BH correction across all tests."""
    frames = []
    for d in datasets:
        res = pair_stats(classify(d["annot"], classes), int(d["n"]))
        res.insert(0, "source", d["source"])
        res.insert(0, "species", d["species"])
        frames.append(res)
    res = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    if res.empty:
        return res
    res["q"] = false_discovery_control(res["p"].to_numpy(), method="bh")
    res["call"] = np.where(res["q"] < alpha, res["direction"], "n.s.")
    low = (res["direction"] == "avoidance") & (res["E"] < min_e)
    res.loc[low, "call"] = f"undetectable (E<{min_e:g})"
    return res


def demo_datasets(seed=0):
    """Synthetic data with planted signals: ring avoids znf, krab co-occurs with znf, bromo with phd."""
    rng = np.random.default_rng(seed)
    names = dict(znf="zf-C2H2 Zinc finger, C2H2 type", ring="zf-RING_2 RING finger domain",
                 krab="KRAB KRAB box", bromo="Bromodomain", phd="PHD PHD-finger", wd40="WD40 WD domain",
                 tpr="TPR_1 Tetratricopeptide repeat", hmg="HMG_box HMG (high mobility group) box",
                 bzip="bZIP_1 bZIP transcription factor", wing="Winged helix DNA-binding domain",
                 brct="BRCT BRCA1 C Terminus domain")
    sizes = dict(znf=800, ring=300, krab=150, bromo=60, phd=120, wd40=400, tpr=200, hmg=50,
                 bzip=60, wing=100, brct=40)
    out = []
    for species, n, scale in [("Demo metazoan", 20000, 1.0), ("Demo protist", 5000, 0.15)]:
        idx = {c: set(rng.choice(n, max(1, int(s * scale)), replace=False).tolist()) for c, s in sizes.items()}
        idx["ring"] -= idx["znf"]
        idx["krab"] = set(rng.choice(sorted(idx["znf"]), len(idx["krab"]), replace=False).tolist())
        idx["phd"] |= set(sorted(idx["bromo"])[: len(idx["bromo"]) // 2])
        rows = [(f"{species.split()[1]}_{g}", names[c]) for c, genes in idx.items() for g in genes]
        out.append(dict(species=species, source="synthetic", n=n,
                        annot=pd.DataFrame(rows, columns=["gene_id", "domain"])))
    return out
