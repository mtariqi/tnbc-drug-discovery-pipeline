# data_loaders/

## Real code (not placeholders)

| File | Status |
|---|---|
| `kinase_data_fetchers.py` | DGIdb/FAERS fetchers. **Verified**: all 3 offline smoke tests pass (mocked network responses) — real network calls not exercised in this sandbox (no network access here). |
| `string_network_builder.py` | STRING API + Louvain communities. **Verified**: all 3 smoke tests pass, including a real correctness check (two tight triangles correctly split into different communities, isolated node gets its own, weak bridge edge doesn't falsely merge them). |
| `tcga_brca_survival_pipeline.py` | Cox/log-rank/RMST. **Verified**: the no-external-dependency logic (survival-time derivation, median stratification, manual log-rank fallback) passes real correctness checks — a null case gives an unremarkable p-value (0.76) and a deliberately planted effect is correctly detected (p≈0). `lifelines` isn't installed in this sandbox, so Cox regression, the Schoenfeld PH test, and RMST specifically were not exercised — install `lifelines` to check those. |
| `gdc_local_data_loader.py` | **Known gap, not silently worked around** — its own docstring says it was tested against a synthetic GDC directory tree, but the code that builds that tree isn't in this file. Crashes with `StopIteration` if run as-is. See `docs/limitations.md` for what to do instead of guessing at the real GDC file format. |

## Real confirmed data path

`~/rtk_nrtk_tnbc/data/` (raw/ and processed/ subfolders) — see `data/README.md` for the exact
files under it, confirmed from actually running `validate_pipeline_inputs.py` and
`run_discovery.py` against it.
