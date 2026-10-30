#!/usr/bin/env python3
"""
Module: app
Purpose: Streamlit UI for genome-wide TF protein domain co-occurrence (Parikesit, Prohaska & Stadler 2011).
Author: Dr. Arli Aditya Parikesit
Date: 2026
Parameters:
    - AUGUSTUS_BIN (env): AUGUSTUS executable, default 'augustus' on PATH
References:
    - Parikesit AA, Prohaska SJ, Stadler PF (2011). Leipzig Research Festival for Life Science.
"""
import tempfile
from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st

from cooc.core import DEFAULT_CLASSES, analyze, count_fasta, demo_datasets, read_annotation
from cooc.pipeline import GP_LOCK, augustus_path, augustus_proteins, hmm_annotate, run_augustus

st.set_page_config(page_title="TF Domain Co-occurrence", layout="wide")
st.title("Genome-wide protein domain co-occurrence of transcription factors")
st.caption("E(x,y) = X·Y/n versus observed F, Poisson test, after Parikesit, Prohaska & Stadler (2011).")

ss = st.session_state
MAX_GENOME_BYTES = 100 * 1024 * 1024
ss.setdefault("datasets", [])

with st.sidebar:
    st.header("Settings")
    alpha = st.number_input("alpha (on BH q-value)", 0.0001, 0.5, 0.05, format="%.4f")
    min_e = st.number_input("Minimum E for detectable avoidance", 0.0, 100.0, 1.0)
    st.subheader("Domain classes (first match wins)")
    classes_df = st.data_editor(pd.DataFrame(DEFAULT_CLASSES, columns=["class", "regex"]),
                                num_rows="dynamic", width="stretch")
    if st.button("Load synthetic demo"):
        ss.datasets = demo_datasets()
    if st.button("Clear datasets"):
        ss.datasets = []

tab_data, tab_gp, tab_pairs, tab_cmp, tab_method = st.tabs(
    ["Data", "Gene prediction (AUGUSTUS)", "Pair table", "Species comparison", "Method"])

with tab_data:
    st.markdown("Upload a domain annotation per species: InterProScan TSV, hmmscan `--domtblout`, "
                "or a two-column `gene_id,domain` table. **n must include genes without domains**; "
                "an underestimated n inflates E and biases calls toward avoidance.")
    with st.form("upload", clear_on_submit=True):
        up = st.file_uploader("Annotation file", type=["tsv", "csv", "txt", "domtbl", "domtblout"])
        fasta = st.file_uploader("Optional proteome FASTA (sets n)", type=["fa", "fasta", "faa"])
        c1, c2, c3 = st.columns(3)
        species = c1.text_input("Species", max_chars=100)
        source = c2.selectbox("Source", ["GP", "SF", "other"])
        n_in = c3.number_input("Total genes n (0 = from FASTA)", 0, 10_000_000, 0)
        if st.form_submit_button("Add dataset") and up and species:
            try:
                annot = read_annotation(up.getvalue(), up.name)
                n = n_in or (count_fasta(fasta.getvalue()) if fasta else 0)
                if not n:
                    raise ValueError("Give n or a proteome FASTA")
                ss.datasets.append(dict(species=species, source=source, n=n, annot=annot))
            except Exception as e:
                st.error(f"Could not add dataset: {e}")
    if ss.datasets:
        st.dataframe(pd.DataFrame([dict(species=d["species"], source=d["source"], n=d["n"],
                                        hits=len(d["annot"]), genes=d["annot"]["gene_id"].nunique())
                                   for d in ss.datasets]), width="stretch")

with tab_gp:
    if augustus_path() is None:
        st.warning("AUGUSTUS not found. Install it and set the AUGUSTUS_BIN (and AUGUSTUS_CONFIG_PATH) "
                   "environment variables before starting the app.")
    else:
        st.markdown("Genome FASTA -> AUGUSTUS gene models -> proteins -> pyhmmer `hmmsearch` with your HMM "
                    "library (e.g. a Pfam-A subset). Fused gene models can create false co-occurrences.")
        with st.form("gp"):
            genome = st.file_uploader("Genome FASTA", type=["fa", "fasta", "fna"])
            hmm = st.file_uploader("HMM library (.hmm, HMMER3)", type=["hmm"])
            c1, c2, c3 = st.columns(3)
            aug_species = c1.text_input("AUGUSTUS --species", "human", max_chars=64)
            label = c2.text_input("Species label", max_chars=100)
            evalue = c3.number_input("Inclusion E-value", 1e-30, 1.0, 1e-3, format="%.1e")
            gathering = st.checkbox("Use gathering (GA) cutoffs instead of E-value")
            run = st.form_submit_button("Run")
        if run and genome and hmm and label:
            if not GP_LOCK.acquire(blocking=False):
                st.error("Another gene prediction run is in progress; try again later.")
                st.stop()
            try:
                if genome.size > MAX_GENOME_BYTES or not genome.getvalue().lstrip().startswith(b">"):
                    raise ValueError("Genome must be a FASTA file of at most 100 MB")
                with tempfile.TemporaryDirectory() as tmp, st.spinner("Running AUGUSTUS and hmmsearch..."):
                    g, h = Path(tmp, "genome.fa"), Path(tmp, "lib.hmm")
                    g.write_bytes(genome.getvalue())
                    h.write_bytes(hmm.getvalue())
                    gff = run_augustus(g, aug_species)
                    prots = augustus_proteins(gff)
                    if not prots:
                        raise ValueError("AUGUSTUS predicted no proteins")
                    annot = hmm_annotate(prots, h, evalue=evalue, gathering=gathering)
                ss.datasets.append(dict(species=label, source="GP", n=len(prots), annot=annot))
                ss.last_gff = gff
                st.success(f"{len(prots)} proteins predicted, {annot['gene_id'].nunique()} with domain hits.")
            except Exception as e:
                st.error(f"Pipeline failed: {e}")
            finally:
                GP_LOCK.release()
        if "last_gff" in ss:
            st.download_button("Download AUGUSTUS GFF3", ss.last_gff, "augustus.gff3")

res = pd.DataFrame()
if ss.datasets:
    try:
        res = analyze(ss.datasets, classes_df[["class", "regex"]].dropna().itertuples(index=False),
                      alpha=alpha, min_e=min_e)
    except ValueError as e:
        st.error(str(e))

with tab_pairs:
    if res.empty:
        st.info("Add a dataset or load the synthetic demo.")
    else:
        calls = st.multiselect("Show calls", sorted(res["call"].unique()), sorted(res["call"].unique()))
        st.dataframe(res[res["call"].isin(calls)], width="stretch")
        st.download_button("Download CSV", res.to_csv(index=False), "cooccurrence.csv", "text/csv")

with tab_cmp:
    if not res.empty:
        res["dataset"] = res["species"] + " (" + res["source"] + ")"
        mat = res.pivot_table(index="pair", columns="dataset", values="log2_F_over_E")
        sig = res.assign(m=res["call"].map(lambda c: "*" if c in ("avoidance", "co-occurrence") else ""))
        sig = sig.pivot_table(index="pair", columns="dataset", values="m", aggfunc="first").reindex_like(mat)
        lim = max(1.0, float(mat.abs().max().max()))
        fig = px.imshow(mat, color_continuous_scale="RdBu_r", zmin=-lim, zmax=lim, aspect="auto",
                        labels=dict(color="log2((F+0.5)/(E+0.5))"))
        fig.update_traces(text=sig.fillna("").values, texttemplate="%{text}")
        fig.update_layout(height=max(400, 22 * len(mat)))
        st.plotly_chart(fig, width="stretch")
        st.caption("Red: more co-occurrence than expected; blue: avoidance. * = significant (BH q < alpha).")

with tab_method:
    st.markdown(r"""
For each pair of domain classes $x, y$ in a genome with $n$ genes:

- $X$, $Y$: genes containing a domain of class $x$ or $y$; $F$: genes containing both.
- Expected co-occurrence $E(x,y) = XY/n$.
- $F < E$: avoidance, $P = \Pr[\mathrm{Pois}(E) \le F]$. $F \ge E$: co-occurrence, $P = \Pr[\mathrm{Pois}(E) \ge F]$.
- Benjamini-Hochberg correction over all pairs and datasets. Effect size $\log_2((F+0.5)/(E+0.5))$.
- When $E < 1$, avoidance cannot be significant and is reported as undetectable (small genomes).

Each domain hit is assigned to the first class whose regex matches, so a RING or PHD finger is not also
counted as a generic zinc finger. Edit the class table to match your HMM library.

**Deviations from the 2011 paper:** AUGUSTUS replaces GENSCAN; any HMMER3 library (e.g. Pfam) replaces
SUPERFAMILY SCOP HMMs; BH correction added. Results are computational hypotheses that need experimental
validation (e.g. cDNA/RT-PCR confirmation of predicted fused gene models, co-immunoprecipitation).
""")
