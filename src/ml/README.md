# ml/

## Built and verified: dependency_ml_model.py

DepMap CRISPR-dependency prediction model — the real Phase 2 "learned model" this project's
grant proposal roadmap called for. Grouped cross-validation by cell line (not random k-fold),
gradient-boosted trees, plus the gene-identity-baseline diagnostic already described in
`results/reports/Heuristic_vs_ML_Report.docx`.

**Verified in this session, two ways:**
1. Its own smoke test (planted signal, `build_feature_matrix()` + `train_and_evaluate()`):
   passes cleanly (R²=0.83, Spearman ρ=0.89, correctly identifies the two truly-informative
   features as most important).
2. A second, harder test I built specifically for the gene-identity-baseline functions
   (`compute_gene_identity_baseline()`, `train_and_evaluate_with_gene_baseline()`), which
   weren't covered by the file's own smoke test: designed a synthetic scenario mirroring the
   real DepMap finding (strong gene-identity baseline, weak per-sample signal) and confirmed
   the code reproduces the same qualitative pattern reported in `Heuristic_vs_ML_Report.docx`
   — gene-identity baseline (R²=0.98) far outperforms the weak multi-omic-only model
   (R²=-0.27), and combining both barely beats gene-identity alone.

## Next step

Run this against the real 25-cell-line DepMap feature matrix (needs `depmap_multiomic_loader.py`
and `depmap_supplemental_loader.py`, referenced in this file's docstrings but not yet part of
this repo — check `tnbc_ml_trained_pipeline`'s original repo for these) to reproduce the real
R²=0.337 (gene-identity) vs. R²=0.339 (combined) finding directly, rather than only on
synthetic data.
