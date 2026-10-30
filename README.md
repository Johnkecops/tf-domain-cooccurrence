# tf-domain-cooccurrence

A Streamlit app that measures genome-wide protein domain co-occurrence and avoidance in transcription factors. It follows Parikesit, Prohaska & Stadler (2011), *Quantitative Measurement of Genome-wide Protein Domain Co-occurrence of Transcription Factors*.

For each pair of domain classes (znf, bZIP, winged helix, bromo, BRCT, KRAB, HMG box, WD40, PHD, RING, TPR) the app computes the expected co-occurrence E = X·Y/n, compares it with the observed count F, and tests the deviation with a one-sided Poisson test. P-values are corrected with Benjamini-Hochberg across all pairs and datasets. A heatmap of log2(F/E) compares species and data sources (gene prediction "GP" vs. database annotation "SF"), in the spirit of the paper's figures.

## Inputs

- Domain annotation per species: InterProScan TSV, HMMER `hmmscan --domtblout`, or a two-column `gene_id,domain` table, plus the total gene count n (or a proteome FASTA to count it). n has to include genes without domains. If it is too small, E is inflated and calls lean toward avoidance.
- UniProt reference proteome (real data, fetched on demand): enter a proteome ID such as UP000005640 (human). The app downloads the proteome's gene2acc file and keeps one canonical protein per gene, takes the UniProtKB Pfam cross-references of those proteins, and names them with Pfam-A.clans.tsv. n is the number of genes in gene2acc, so genes without Pfam domains are counted. Proteins whose gene is unknown are left out of both n and the annotation. This is database annotation (the "SF" side of the paper), and it needs network access to rest.uniprot.org, ftp.uniprot.org and ftp.ebi.ac.uk. A large proteome takes a few minutes, because UniProt streams can be slow. Results are cached per server process.
- Genome FASTA (optional): AUGUSTUS predicts the genes and pyhmmer `hmmsearch` annotates the predicted proteins against an HMMER3 library you supply, such as a Pfam-A subset. This replaces the GENSCAN + SUPERFAMILY setup of the original work.
- A synthetic demo with planted signals (krab-znf and bromo-phd co-occurrence, ring-znf avoidance). It is not real biology.

## Run

```bash
pip install -r environments/requirements.txt
# gene prediction tab: AUGUSTUS is found on PATH or at the Homebrew location
# (/opt/homebrew/opt/augustus/bin/augustus, tested with 3.5.0). Override with:
# export AUGUSTUS_BIN=/path/to/augustus AUGUSTUS_CONFIG_PATH=/path/to/augustus/config
streamlit run app.py
python -m pytest -q tests
```

Tested with Python 3.11.5, streamlit 1.57.0, pandas 2.3.3, scipy 1.17.1, plotly 6.7.0, pyhmmer 0.12.3, requests 2.31.0.

## Caveats

- If E < 1 (typical for small genomes), avoidance cannot reach significance, so the app reports it as undetectable.
- A domain hit goes to the first class whose regex matches, which keeps RING/PHD/KRAB from also being counted as generic znf. Adjust the class table in the sidebar to fit your HMM library.
- De novo gene models can fuse neighbouring genes and produce false co-occurrences. Treat every call as a hypothesis and confirm it experimentally, for example with cDNA/RT-PCR of the predicted gene structure.

## References

- Parikesit AA, Prohaska SJ, Stadler PF (2011). Quantitative Measurement of Genome-wide Protein Domain Co-occurrence of Transcription Factors. Leipzig Research Festival for Life Science, University of Leipzig.
- Parikesit AA, Stadler PF, Prohaska SJ (2010). Quantitative comparison of genomic-wide protein domain distributions. German Conference on Bioinformatics, LNI P-173:93-102.
- Stanke M, Morgenstern B (2005). AUGUSTUS: a web server for gene prediction in eukaryotes that allows user-defined constraints. Nucleic Acids Res 33:W465-W467. doi:10.1093/nar/gki458
- The UniProt Consortium (2025). UniProt: the Universal Protein Knowledgebase in 2025. Nucleic Acids Res 53:D609-D617. doi:10.1093/nar/gkae1010
- Blum M et al. (2025). InterPro: the protein sequence classification resource in 2025. Nucleic Acids Res 53:D444-D456. doi:10.1093/nar/gkae1082
- Larralde M, Zeller G (2023). PyHMMER: a Python library binding to HMMER for efficient sequence analysis. Bioinformatics 39(5):btad214. doi:10.1093/bioinformatics/btad214

## Security

The /security-audit harness has been run on this codebase, and the vulnerabilities it found were patched (regex timeout against ReDoS, localhost-only binding, upload and genome size caps, a single-run lock for AUGUSTUS, AUGUSTUS stderr kept in the server log).
The app has no authentication, so a shared deployment needs an authenticating reverse proxy in front of it.

## License

GPLv3, see [LICENSE.md](LICENSE.md).

AI Assistance Disclaimer: This codebase was developed with the assistance of Claude Code. While the AI provided code generation, debugging, and structural support, the human developer maintains full responsibility for reviewing, testing, and maintaining all content and functionality
