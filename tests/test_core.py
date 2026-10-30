import pandas as pd
import pytest

from cooc.core import DEFAULT_CLASSES, analyze, classify, demo_datasets, pair_stats, read_annotation
from cooc.pipeline import augustus_proteins
from cooc import uniprot
from cooc.uniprot import canonical_genes, pfam_annotation, pfam_names


def test_expectation_and_direction():
    # 10 znf genes, 10 ring genes, overlap 0, n = 100 -> E = 1, avoidance
    gc = pd.DataFrame([(f"g{i}", "znf") for i in range(10)] + [(f"g{i}", "ring") for i in range(10, 20)],
                      columns=["gene_id", "cls"])
    r = pair_stats(gc, 100).iloc[0]
    assert r["E"] == pytest.approx(1.0) and r["F"] == 0 and r["direction"] == "avoidance"
    assert r["p"] == pytest.approx(0.36787944, rel=1e-6)  # Pois(0; 1) = e^-1
    with pytest.raises(ValueError):
        pair_stats(gc, 5)


def test_first_match_wins():
    annot = pd.DataFrame([("a", "zf-RING_2 RING finger"), ("b", "zf-C2H2"), ("c", "HMG-CoA reductase")],
                         columns=["gene_id", "domain"])
    got = dict(classify(annot, DEFAULT_CLASSES).values)
    assert got == {"a": "ring", "b": "znf"}


def test_demo_recovers_planted_signals():
    res = analyze(demo_datasets()).set_index(["species", "pair"])
    assert res.loc[("Demo metazoan", "krab-znf"), "call"] == "co-occurrence"
    assert res.loc[("Demo metazoan", "ring-znf"), "call"] == "avoidance"
    assert res.loc[("Demo metazoan", "bromo-phd"), "call"] == "co-occurrence"


def test_parsers():
    ipr = "P1\tmd5\t300\tPfam\tPF00096\tZinc finger, C2H2 type\t1\t20\t1e-5\tT\tdate\tIPR013087\tZnf C2H2\n"
    assert read_annotation(ipr).values.tolist() == [["P1", "PF00096 Zinc finger, C2H2 type Znf C2H2"]]
    assert read_annotation("gene_id,domain\ng1,KRAB\n").values.tolist() == [["g1", "KRAB"]]
    gff = ("chr1\tAUGUSTUS\tmRNA\t1\t90\t.\t+\t.\tID=g1.t1;Parent=g1\n"
           "# protein sequence = [MKV\n# LLA]\n")
    assert augustus_proteins(gff) == {"g1.t1": "MKVLLA"}


def test_redos_pattern_times_out():
    annot = pd.DataFrame([("g1", "a" * 40 + "!")], columns=["gene_id", "domain"])
    with pytest.raises(ValueError, match="timed out"):
        classify(annot, [("evil", r"(a|aa)+$")])


def test_uniprot_parsers():
    g2a = "-\tQ0\t-\nG1\tP1\tG1\nG1\tA1\tG1\nRT1\tR1\tG1\nG2\tP2\tG2\n"
    canon = canonical_genes(g2a)
    assert canon == {"P1": "G1", "P2": "G2"}  # A1, read-through R1, unknown-gene Q0 dropped
    names = pfam_names("PF00096\tCL0361\tC2H2-zf\tzf-C2H2\tZinc finger, C2H2 type\n")
    tsv = "Entry\tPfam\nP1\tPF00096;PF01352;\nA1\tPF00096;\nQ0\tPF00096;\nP2\t\n"
    got = pfam_annotation(tsv, canon, names).values.tolist()
    assert got == [["G1", "zf-C2H2 PF00096 Zinc finger, C2H2 type"], ["G1", "PF01352"]]


def test_truncated_gzip_is_retried(monkeypatch):
    import gzip
    full = gzip.compress(b"Entry\tPfam\nP1\tPF00096;\n")
    bodies = [full[:-8], full]  # first response cut off before the gzip trailer

    class Resp:
        def __init__(self, content):
            self.content = content

        def raise_for_status(self):
            pass

    monkeypatch.setattr(uniprot.requests, "get", lambda url, timeout: Resp(bodies.pop(0)))
    monkeypatch.setattr(uniprot, "RETRY_WAIT", 0)
    assert uniprot._get("https://example.org/x") == "Entry\tPfam\nP1\tPF00096;\n"
    bodies[:] = [full[:-8]] * 3
    with pytest.raises(RuntimeError, match="after 3 attempts"):
        uniprot._get("https://example.org/x")
