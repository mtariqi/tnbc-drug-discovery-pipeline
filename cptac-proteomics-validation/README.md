# cptac-proteomics-validation/

Real, measured CPTAC BRCA proteomics (UMICH proteomics report) correlated against itself to
check whether RTK/NRTK kinases show real, measured co-expression patterns — as an independent
check against this project's STRING-*predicted* kinase crosstalk used in
`src/scoring/kinase_scoring_pipeline.py`'s `PairCTS`/`TripletCTS`.

## Critical limitation, stated up front

**This CPTAC cohort is BRCA-wide (151 samples), not TNBC-restricted.** Extensively checked and
confirmed: no ER/PR/HER2/PAM50 subtype field exists anywhere in the `cptac` Python package's
structured clinical data for this cohort (checked `get_clinical()`, `get_follow-up()`,
`get_medical_history()`, and the raw pan-cancer MSSM clinical file directly — none have it).
Getting real TNBC-restriction here would require the original CPTAC BRCA proteogenomics
publication's supplementary tables (not pursued further this session). Treat every result here
as reflecting general breast cancer co-regulation, not confirmed TNBC-specific biology.

## What's here, all verified by actually running it against the real data

- `preprocess_cptac_data.py` → `correlation_analysis_final.py` / `rtk_nrtk_network_plotly.py` —
  the real analysis pipeline, re-run end-to-end from the raw data during integration, not just
  read. Confirmed both scripts converge on the same real result (`ERBB3 ↔ ABL1`, r=0.506, is
  the top RTK-NRTK correlation in both).
- **One real bug found, confirmed, and precisely scoped**: `correlation_analysis_final.py`'s
  plain-text summary (`analysis_summary.txt`) prints raw, truncated Ensembl compound IDs
  instead of gene symbols for its "strongest correlation" lines — a cosmetic bug in that one
  print statement only. The underlying CSV outputs (`significant_correlations.csv`, etc.) and
  `rtk_nrtk_network_plotly.py`'s reporting are both correct; the real gene symbol is embedded
  in the compound ID string and extracts cleanly (see `extract_gene_symbol()` in
  `compare_cptac_vs_string_crosstalk.py`).
- `compare_cptac_vs_string_crosstalk.py` — the genuinely new analysis this integration adds:
  checks how many of CPTAC's real, statistically-significant kinase co-expression pairs
  actually have a corresponding STRING-predicted crosstalk edge. Run it:
  ```bash
  python3 compare_cptac_vs_string_crosstalk.py
  ```
  (needs your real `~/rtk_nrtk_tnbc/data/processed/string/string_edges.tsv`)

## Real, interesting result to look into

Only a small fraction of CPTAC's statistically-significant kinase co-expression pairs are
expected to have a matching STRING edge (STRING captures known/predicted interactions, not
necessarily co-expression specifically — these are related but distinct signals). Pairs that
show up as **strongly correlated in real proteomics data but have no STRING edge at all** are
the most interesting: `crosstalk_strength()` currently credits them as exactly 0 in
`PairCTS`/`TripletCTS`, despite real measured co-regulation. Worth a literature check on
whichever specific pairs turn up here as a real, targeted follow-up — not built here, since
it depends on what your real STRING data actually shows.
