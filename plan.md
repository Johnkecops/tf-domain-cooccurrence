# Plan: TF Domain Co-occurrence Streamlit App

Source paper: Parikesit, Prohaska & Stadler (2011), *Quantitative Measurement of Genome-wide Protein Domain Co-occurrence of Transcription Factors*, Leipzig Research Festival for Life Science.

## 1. What the paper does (to reproduce)

- Domains of interest. Seven TF DNA-binding classes: zinc finger (znf), leucine zipper (bzip), winged helix (wing), bromo, BRCT, KRAB, HMG box. Non-DNA-binding partners of znf: WD40, PHD, RING, TPR.
- Two data sources per species: de novo gene prediction + HMM domain annotation (GP), and SUPERFAMILY database annotation (SF).
- Statistic per domain pair (x, y):
  - X = genes with domain x, Y = genes with domain y, n = total genes
  - Expected co-occurrence E(x, y) = X * Y / n
  - Observed F = genes carrying both x and y
  - F < E means avoidance; F > E means preferred co-occurrence
  - Significance: Poisson with mean E (P < 0.05)
- Caveat from the paper: when E < 1 (small genomes), avoidance cannot be detected.
- 18 eukaryotes (Giardia to Homo sapiens) compared side by side.

## 2. Deviations from the paper (deliberate)

| Paper | App | Reason |
|---|---|---|
| GENSCAN gene prediction | AUGUSTUS (local binary, `--species` user-selected) | User instruction; AUGUSTUS is maintained and open source |
| SCOP HMMs (SUPERFAMILY) | Any HMM library the user supplies (Pfam-A subset recommended), scanned with pyhmmer 0.12.3 | SUPERFAMILY HMMs need a licence download; pyhmmer is already installed, no HMMER binary needed |
| Single Poisson P | One-sided Poisson P in the observed direction + Benjamini-Hochberg FDR across all tests, plus log2(F/E) effect size | Many pairs x species are tested at once |

## 3. Inputs

1. **Domain annotation table (main path, fast).** Upload one or more files, each tagged with species and source (GP / SF / other). Supported formats:
   - InterProScan TSV (column 1 protein ID, columns 5-6 signature accession/description)
   - HMMER `--domtblout`
   - Generic 2-column CSV/TSV `gene_id,domain`
   - Total gene count n is entered or taken from a protein FASTA (proteins without domains still count toward n).
2. **Genome FASTA (full GP path).** AUGUSTUS predicts genes, proteins are extracted from the GFF3 (`--protein=on`), pyhmmer `hmmscan` against an uploaded `.hmm` file produces the domain table, which then feeds path 1.
3. **Demo dataset.** Small synthetic table bundled in `data/demo/`, clearly labelled as synthetic, so the app runs without any downloads.

## 4. Domain classification

Regex keyword map from domain name/description to class (e.g. `zf-|zinc finger` -> znf, `KRAB` -> krab, `WD40|WD repeat` -> wd40). Editable in the sidebar as a table, so users can tune it for their HMM library. Genes are counted once per class (presence/absence), as in the paper.

## 5. App layout (single `app.py`)

- Sidebar: data input, domain class map, alpha, FDR toggle, minimum E threshold.
- Tab "Pair table": X, Y, n, E, F, log2(F/E), P, q, call (avoidance / co-occurrence / n.s. / undetectable E<1), CSV download.
- Tab "Species comparison": heatmap of log2(F/E) for domain pair x species, significance marked, GP and SF side by side (mirrors paper Figures 3-9).
- Tab "Gene prediction": AUGUSTUS + pyhmmer runner with log and downloads.
- Tab "Method": formulas, citation, caveats (predictions are hypotheses; GP artefacts like fused gene models inflate co-occurrence).

## 6. File layout

```
app.py                 Streamlit UI
cooc/core.py           parsing, classification, E/F/Poisson/BH (pure functions)
cooc/pipeline.py       AUGUSTUS + pyhmmer wrapper (subprocess, no shell)
data/demo/             synthetic demo table
tests/test_core.py     one small test file: E formula, Poisson direction, BH, parsers
environments/requirements.txt
README.md, LICENSE.md (GPLv3), plan.md, .gitignore
```

## 7. Security considerations (input to the audit)

- `subprocess.run([...])` with argument lists only, never `shell=True`; species name validated against `augustus --species=help` / regex `^[A-Za-z0-9_]+$`.
- Uploads written to a per-session temp dir, size-capped (Streamlit `maxUploadSize`), deleted after use.
- AUGUSTUS runtime timeout; FASTA validated before execution.
- No network calls at runtime.

## 8. Post-build steps (from instruction.md)

1. Run `/security-audit` on the codebase and patch findings.
2. Add two lines to README.md stating the security-audit harness was run and vulnerabilities patched.
3. Write LICENSE.md (GPLv3).
4. `git init`, commit, push to a new repository under https://github.com/Johnkecops.

## 9. Open items

- AUGUSTUS binary was **not found** on PATH, Homebrew, or conda envs at planning time. The app takes the binary path as a setting (`AUGUSTUS_BIN` env var or sidebar field) and disables the GP tab with a clear message if absent.
- Results from GP mode are computational hypotheses; experimental validation (e.g. RT-PCR of predicted fusion genes, co-IP) is required before biological claims.

## Reference

- Parikesit AA, Prohaska SJ, Stadler PF (2011). Quantitative Measurement of Genome-wide Protein Domain Co-occurrence of Transcription Factors. Leipzig Research Festival for Life Science, University of Leipzig, 17 Dec 2010 (poster abstract).
- Parikesit AA, Stadler PF, Prohaska SJ (2010). Quantitative comparison of genomic-wide protein domain distributions. German Conference on Bioinformatics, LNI P-173:93-102.
- Stanke M, Morgenstern B (2005). AUGUSTUS: a web server for gene prediction in eukaryotes. Nucleic Acids Res 33:W465-W467. doi:10.1093/nar/gki458
- Larralde M, Zeller G (2023). PyHMMER: a Python library binding to HMMER. Bioinformatics 39(5):btad214. doi:10.1093/bioinformatics/btad214
