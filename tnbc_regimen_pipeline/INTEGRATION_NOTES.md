# Integration Notes — tnbc_regimen_pipeline merge

## What this is, and how it relates to src/regimen/ and src/discovery/

This is a **real, separately-versioned (v3.0.0), installable Python package** — pyproject.toml,
proper package structure, its own test suite — that appears to be a matured, productionized
version of the ad-hoc scripts in `src/regimen/` and `src/discovery/` (`cohort_wide_regimen_analysis.py`
and `agentic_regimen_discovery.py` specifically). It adds a CLI, YAML-driven job configs,
before/after discovery diffing, provenance stamping, caching, and parallelized discovery —
none of which the ad-hoc scripts had.

**Design discipline it follows, worth preserving:** it deliberately does NOT reimplement
MDCOE/HCOS scoring — `resolve_tp53_fn`, `DrugGraph`, `SynergyNet`, `hcos_fn`, `mdcoe_fn` are
taken as plugin callables (a `--plugin module:function` CLI arg, or `plugin:` in a YAML job
config), keeping the real, separately-validated `mdcoe.py` as the single source of truth for
scoring logic rather than duplicating it here.

**Not yet resolved:** `src/regimen/cohort_wide_regimen_analysis.py` and
`src/discovery/agentic_regimen_discovery.py` (the older, ad-hoc versions) still exist
side-by-side with this package's `cohort/cohort_wide_regimen_analysis.py` and
`discovery/agentic_regimen_discovery.py`. Decide whether to keep both (the old ones as a
lighter-weight/scriptable alternative) or retire the older versions now that this package
supersedes them — don't just let both silently drift.

## Verified in this session, not taken on the README's word

1. **Test suite: 21/21, independently re-confirmed.** The package's own README states it was
   verified via a manual harness (no pytest/network in the original build sandbox) and asks
   you to run real `pytest tests/` to confirm. Ran it with my own already-fixed harness
   (`_manual_test_harness.py` — the same one that caught and fixed the conftest/class-based-test
   gaps in `tnbc-genomics-agent`'s suite) instead of the original, unfixed harness this package
   shipped with: **21 passed, 0 failed**, matching the README's claim exactly. (The original
   flat-function test style here didn't actually need the class-based-test fix — but using the
   fixed harness anyway costs nothing and is more defensible.)

2. **Full CLI pipeline, run genuinely end-to-end** — not just the internal test suite's mocks.
   Built a real gzipped MAF file (2 patients, one being the actual TCGA-AO-A128 gene set:
   EGFR/PTEN/FLT1/TP53/ERBB2/TYK2), a real plugin module importing the actual verified
   `mdcoe.py`, and ran `tnbc-pipeline full-pipeline` for real:
   - **Reproduced the headline result to floating-point precision**: `TCGA-AO-A128` →
     `afatinib + alpelisib + trastuzumab`, HCOS = 0.45000000000000007.
   - **Confirmed real network-failure handling, not just a mocked success path**: the
     discovery step made an actual live HTTPS call to PubMed's E-utilities (blocked here —
     no real network access in this sandbox — got a genuine 403), and the pipeline logged the
     failure, excluded that one gene from discovery results, and completed successfully rather
     than crashing. This is a meaningfully different (and better) check than only running
     against the test suite's mocked network responses.
   - Output artifacts (`cohort_result.csv`, `regimen_frequency.csv`, `merged_gene_drugs.json`,
     `provenance.json`, `regimen_diff.csv`, `discovery_result.csv`) all generated correctly,
     each stamped with `_pipeline_version`, `_git_commit`, `_run_timestamp_utc` for real
     reproducibility tracking.

## Known open items (stated in the package's own README, carried forward here)

- `resolve_tp53_fn`'s default (`Nonsense_Mutation` for every patient) is only confirmed for
  the one previously-inspected patient (TCGA-AO-A128) — verify each real cohort patient's
  actual TP53 classification before trusting cohort-wide TP53-arm results at scale.
- Schema validation is hand-written, not `pydantic` (unavailable when built) — fine for the
  current shallow dict shape, worth revisiting if the schema grows.
- The gene→drug→regimen flow SVG is hand-rolled, not a true Sankey (`plotly` unavailable).
- `LICENSE` is MIT as a placeholder — confirm this is actually right for your institution
  before publishing.

## Real, confirmed data-path note carried over from assemble_tcga_brca_data.R

The R data-assembly script (now in `scripts/r/assemble_tcga_brca_data.R`) contains its own
self-flagged discrepancy: `data_dir <- "~/rtk_nrtk_tnbc/data/raw/tcga_brca/GDCdata"`, with an
inline comment noting a screenshot showed `.../raw/tcga_brca/...` while some other script used
`.../raw/tcga/...`. **Not resolved here** — confirm which is the real path on your machine
before running this script, rather than assume the comment's already-corrected version is
necessarily the final one.
