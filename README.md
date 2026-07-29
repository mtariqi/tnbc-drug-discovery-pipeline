![banner](docs/banner.png)

# TNBC Drug Discovery Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![PyTorch](https://img.shields.io/badge/PyTorch-GNN%20models-ee4c2c.svg)](https://pytorch.org/)
[![Status](https://img.shields.io/badge/status-manuscript%20in%20submission-orange.svg)]()
[![DOI](https://img.shields.io/badge/DOI-pending%20Zenodo%20archive-lightgrey.svg)]()
# TNBC Kinase-Redundancy, Drug-Regimen, and Cohort-Scale ML Project

Consolidated repository combining three previously-separate GitHub repos
(`tnbc-kinase-scoring-pipeline`, `TNBC-drug-regimen-discovery`, `tnbc_ml_trained_pipeline`)
plus new cohort-scale (n≈1000, full TCGA-BRCA) ML work, into one structure.

## Project scope

A real-data computational pipeline for RTK/NRTK target prioritization and drug-regimen
discovery in triple-negative breast cancer, extended from a single focal-patient case
(TCGA-AO-A128) to four cohort-scale ML tracks:

| Track | Question | Population | Status |
|---|---|---|---|
| A — Survival/prognosis | Patient-level risk prediction | Full BRCA cohort, n≈1000-1098, subtype as feature | **Built & verified** (`src/survival/patient_level_survival_model.py`) |
| B — Treatment response | Predict pCR/response | — | **Not attempted** — no response label in this data; see `docs/limitations.md` |
| C — TNBC subtyping | Lehmann-style molecular subtypes | TNBC-only, n≈150-200 | **Built & run on real data (1095 patients)** — code verified correct via synthetic ground truth, but real-cohort clustering gives consistently weak separation (silhouette 0.02-0.07) across three independently-tried methodologies. Likely reflects genuine TNBC subtype overlap in bulk RNA-seq, not a fixable bug — see `docs/limitations.md` |
| D — Biomarker/drug-target discovery | Which kinases/regimens generalize across the cohort | TNBC-only, n≈150-200 | **Built & verified end-to-end** — real `mdcoe.py` genuinely reproduces the HCOS=0.450 headline result and wires correctly into `cohort_wide_regimen_analysis.py` (see `src/regimen/README.md`) |

`src/regimen/`, `src/scoring/`, and `src/data_loaders/` now contain **real, uploaded source**
(not placeholders) from your two original repos — see each folder's own README for exactly
what was verified vs. what still needs your real `~/rtk_nrtk_tnbc/data/` files to run for real.


## Four subprojects in this repo

This repo now contains **four independently-built efforts** addressing related but distinct
angles on TNBC kinase-redundancy and combination therapy — nothing here has been
cross-validated against the others except where explicitly noted, and each subproject's own
docs say so:

- **`src/` + `data/`** — the CTS/HCOS/DepMap cohort-scale project (this README's main focus below)
- **`tnbc-genomics-agent/`** — a self-contained, single-patient VCF-to-Claude-narrative pipeline
  (10-kinase curated DB, redundancy/bypass-risk scoring, AI synthesis, companion report in
  `tnbc-genomics-agent/docs/`). Test suite verified passing (31/31).
- **`mc-ore-prototype/`** — a proposed *extension* of `tnbc-genomics-agent` (GNN synergy
  prediction + knowledge graph + multi-agent LLM reasoning), with a runnable Phase-1 scaffold
  notebook. See `mc-ore-prototype/README.md` for a naming-inconsistency note (MC-ORE /
  Hybrid-CORE / CL-MODE across its source documents).
- **`tnbc_regimen_pipeline/`** — a real, installable (v3.0.0) Python package: a matured,
  CLI-driven, provenance-stamped version of `src/regimen/`'s and `src/discovery/`'s ad-hoc
  scripts. **Verified end-to-end in this session**, not just via its own test suite — a real
  CLI run reproduces the headline afatinib+alpelisib+trastuzumab HCOS=0.450 result to
  floating-point precision, and correctly handles a genuine live-network failure (PubMed
  E-utilities 403, no network in this sandbox) by logging and excluding the affected gene
  rather than crashing. See `tnbc_regimen_pipeline/INTEGRATION_NOTES.md`.

`docs/source_repo_readmes/` holds the real READMEs from your two original GitHub repos
(`tnbc-kinase-scoring-pipeline`, `TNBC-drug-regimen-discovery`) — these corrected several
guessed filenames in the `src/scoring/`, `src/data_loaders/`, and (now-removed) `src/utils/`
placeholder READMEs.

- **`cptac-proteomics-validation/`** — real, measured CPTAC BRCA proteomics, checked against
  STRING-predicted kinase crosstalk edges used in `PairCTS`/`TripletCTS`. **Not TNBC-restricted**
  (BRCA-wide, 151 samples) — extensively checked, no subtype metadata exists in this cohort via
  the `cptac` package. See `cptac-proteomics-validation/README.md`.

- **`docs/kinase_panel_audit/`** — a real, cited kinase-panel audit (56 RTK + 32 NRTK = 88,
  vs. the original 54+29=83) plus project background docs. **Found a real, confirmed
  discrepancy**: the actual 90-kinase panel used throughout this entire project includes
  `LMTK2`/`LMTK3`, which this same audit explicitly says to exclude — see `docs/limitations.md`.

## Repository layout

```
tnbc-project/
├── tnbc-genomics-agent/      # separate subproject — see its own README + INTEGRATION_NOTES.md
├── mc-ore-prototype/          # separate subproject — proposed extension of tnbc-genomics-agent
├── tnbc_regimen_pipeline/     # separate subproject — real installable package, verified end-to-end
├── data/
│   ├── raw/            # real, externally-sourced data — see data/README.md for exact
│   │                   # download instructions per source (never hand-edited)
│   └── processed/      # code-generated intermediate files, regenerable from raw/
├── src/
│   ├── data_loaders/   # per-source loaders (TCGA, DepMap, STRING, DGIdb, openFDA, PubMed)
│   ├── scoring/        # CTS / PairCTS / TripletCTS / TripletCTS_v2
│   ├── regimen/        # MDCOE/HCOS beam search, cohort-wide regimen analysis (Track D)
│   ├── survival/       # per-kinase survival test + patient-level risk model (Track A)
│   ├── subtyping/      # TNBC molecular subtyping (Track C)
│   ├── ml/             # DepMap dependency-prediction ML investigation
│   ├── discovery/      # agentic literature-mining pipeline
│   └── utils/          # shared helpers (TP53 resolution, drug graphs, HCOS formula)
├── tests/              # mirrors src/, synthetic-data-with-known-structure smoke tests
├── notebooks/          # exploratory analysis
├── results/
│   ├── figures/
│   ├── tables/
│   └── reports/        # the project's real .docx findings reports
├── docs/                # methodology, data provenance, limitations
└── scripts/             # CLI entry points that wire loaders → src/ modules end-to-end
```

Every subfolder that isn't populated yet has its own `README.md` explaining exactly what
belongs there and which of your existing three repos it should be copied from — nothing here
is guessed at or fabricated.

## Getting started on my own machine

**Using fish shell?** See `docs/FISH_SHELL_GUIDE.md` first — venv activation, env vars
(`ANTHROPIC_API_KEY`, `TNBC_MIN_AF`, etc.), and running each subproject all use different
syntax than bash. `scripts/setup_local_data.fish` is a native fish port of step 1 below.

1. Run `scripts/setup_local_data.sh` (edit the paths at the top first) to copy your already-
   downloaded TCGA/DepMap/STRING/DGIdb files into `data/raw/` in the expected layout.
2. Copy the four not-yet-migrated code groups (see each folder's `README.md`) from your
   existing three repos: `src/scoring/`, `src/ml/`, `src/utils/`, `src/data_loaders/`.
3. `pip install -r requirements.txt`
4. `python tests/../src/survival/patient_level_survival_model.py` to confirm Track A's smoke
   tests still pass in your real environment.
5. Wire real data into `scripts/run_track_a_survival.py` (currently a template with the
   loader calls commented out) and run it.

## Principles this project follows (keep these when extending it)

- Real data only. Synthetic data is confined to smoke tests and never contributes to a
  reported finding.
- State gaps plainly (see `docs/limitations.md`) rather than filling them with an invented
  proxy.
- Verify real file structure (`inspect_clinical_columns()`-style checks) before trusting an
  assumed column name.
- Prefer convergent validation across independent methods over a single scoring approach.
