# regimen/

## Real, verified code (not placeholders)

- `mdcoe.py` — the real MDCOE beam search + HCOS scorer, copied from
  `TNBC-drug-regimen-discovery`. **Verified in this session**: reproduces the headline
  afatinib+alpelisib+trastuzumab, HCOS=0.450 result exactly for the TCGA-AO-A128 gene set.
  Note this is a different/newer version than described in the repo's own README — its
  knowledge base (`GENE_PATHWAYS`, `GENE_DRUGS`, `DRUG_TOXICITY`) is shared with
  `mc-ore-prototype/tnbc_combo_pipeline.ipynb`, and it includes a TP53 zygosity-aware
  `resolve_tp53_drugs()` (mutant vs. deleted vs. NMD-degraded-truncating, not a flat lookup).
- `cohort_wide_regimen_analysis.py` — runs MDCOE across a whole cohort, not just one patient.
  **Verified in this session**: wired end-to-end with the real `mdcoe.py` above (not mocks) —
  patients with identical altered genes correctly converge on the same regimen, different
  gene profiles get different regimens, and gene-less patients get `NaN` rather than a
  fabricated score.
- `my_mdcoe_plugin.py` — dependency-injection glue between the two files above. Path corrected
  to this repo's own layout (was hardcoded to a separate local machine path).
- `build_real_gene_drugs.py` — builds a real, DGIdb-sourced alternative to `mdcoe.py`'s
  hand-curated `GENE_DRUGS`, for convergent-validation testing (see next file).
- `compare_curated_vs_real_gene_drugs.py` — runs the same patient through MDCOE twice (curated
  vs. real-DGIdb candidate lists) and reports both side by side. This is the code behind the
  "convergent validation" claim in `TNBC_Regimen_Recommendation_Final.docx`.
- `run_discovery.py` — real runner script connecting `mdcoe.py`, `build_real_gene_drugs.py`,
  and `src/discovery/agentic_regimen_discovery.py` together. Contains the real, confirmed data
  path convention: `~/rtk_nrtk_tnbc/data/...` (see `data/README.md`).

## Not yet run against real (non-synthetic) data

`kinase_scoring_pipeline.py` (in `src/scoring/`) still only has synthetic demo output verified
— `validate_pipeline_inputs.py` confirms real data isn't present in this sandbox
(`~/rtk_nrtk_tnbc/data/...` files not found here, as expected). Run it on your own machine
where that data actually exists.
