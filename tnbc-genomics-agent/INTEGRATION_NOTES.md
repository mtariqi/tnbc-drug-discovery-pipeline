# Integration Notes — tnbc-genomics-agent merge

## What this subproject is, and how it relates to the rest of the repo

`tnbc-genomics-agent/` is a **separate, self-contained pipeline** from the CTS/HCOS/DepMap
project documented in `results/reports/` and implemented in `src/`. It takes a different
approach to a related question:

| | `src/` (CTS/HCOS project) | `tnbc-genomics-agent/` |
|---|---|---|
| Kinase data | 90-kinase panel, scored from real STRING/DepMap/TCGA-BRCA/DGIdb data | 10-kinase curated dictionary, hand-entered (overexpression frequencies, inhibitors, resistance mechanisms) |
| Patient input | TCGA-BRCA cohort (clinical + expression + MAF) | Single-patient VCF file |
| Scoring | CTS / PairCTS / TripletCTS (data-driven composite scores) | Sigmoid-scaled pathway redundancy + bypass-risk score (formula-based, not fit to external data) |
| Regimen ranking | MDCOE/HCOS beam search over a patient's altered genes | Ranked gene-pair suggestions from `redundancy_analyzer.py`, same general idea, independent implementation |
| Synthesis | Structured reports (this repo's `.docx` files) | Claude API call (`ai_agent.py`) generating a clinician-readable narrative |

**These are not currently cross-validated against each other.** They were built independently
and use different kinase panels (10 vs. 90), different drug dictionaries, and different
scoring formulas. Before treating any output from one as corroborating the other, that
overlap needs to be checked explicitly — it hasn't been yet.

## Discrepancies found while merging, documented rather than silently fixed

1. **README's pytest-free test claim didn't hold as originally written — now genuinely resolved.**
   `README.md` documents `python -m unittest discover tests/ -v` as a no-pytest fallback, but
   `test_vcf_parser.py` and `test_redundancy_analyzer.py` both `import pytest` directly
   (`pytest.raises`, `pytest.fixture`), which fails in a no-network sandbox. A manual
   pytest-emulation harness (`_manual_test_harness.py`) was supplied to solve exactly this, but
   as first written it also failed here: it expected a `tests/conftest.py` (none exists — the
   `analyzer` fixture actually lives inside `test_redundancy_analyzer.py` itself) and only
   discovered flat `test_*` functions, not the class-based `TestRedundancyAnalyzer.test_*`
   methods the real suite actually uses. Both were confirmed by running it, not assumed from
   reading it — it threw `FileNotFoundError` on the missing conftest before any test ran.

   **Fixed and re-verified**: `resolve_and_call()` now searches the test instance, the test
   module itself, and (if present) `conftest.py`, in that order, and `run_test_file()` now
   discovers class-based tests too. Re-run after the fix: **31 tests pass, 0 fail** — one more
   than the 30 the project report states (`test_redundancy_analyzer.py` has 11 tests, not 10;
   a minor count discrepancy in the report, not a test failure). This is now a genuinely
   confirmed result, not a documentation claim taken on faith.

   **Independently re-confirmed with real pytest 9.1.1** (Python 3.13, on the user's own
   machine, not this sandbox): `python -m pytest tests/ -v` → **31 passed, 0 failed**, same
   count, same test names, full agreement with the manual harness. Two independent
   verification paths now agree, not just one tool's self-report.

   Minor unrelated oddity noticed in that run: pytest's `rootdir`/`configfile` resolved to a
   `setup.cfg` in the user's home directory, not this project — harmless here since collection
   and results were unaffected, but worth giving this project its own `pytest.ini` or
   `pyproject.toml` `[tool.pytest.ini_options]` block before adding any real pytest config, so
   it doesn't silently inherit unrelated settings from outside the repo.

2. **`Agentic_AI_Discovery_Findings.docx`-style honesty note on `mdcoe.py`:** the
   `TNBC-drug-regimen-discovery` and `tnbc-kinase-scoring-pipeline` READMEs both list their own
   copy of `mdcoe.py` — the drug-regimen repo's README documents a real beam-search depth bug
   found and fixed there ("bug-fixed & re-verified" badge); the kinase-scoring repo's own
   Known Limitations section separately states its `mdcoe.py` and CTS/TripletCTS scoring "are
   not yet reconciled." **Two copies of the same file, at possibly different bug-fix states,
   is a real drift risk** — worth confirming which copy `my_mdcoe_plugin.py` should actually
   point at (currently hardcoded to the drug-regimen-discovery repo's copy, which per its
   README is the one with the beam-search fix already verified).

2. **The uploaded zip was flat, not in the `modules/` layout its own code expects.**
   `pipeline.py`, `redundancy_analyzer.py`, `vcf_parser.py`, and all three test files import
   from `modules.kinase_db`, `modules.vcf_parser`, etc. — but the zip had every file at the
   top level with no `modules/` folder. This was corrected during the merge (files moved into
   `modules/`), and the full pipeline was re-run end-to-end against `data/sample/patient_1.vcf`
   to confirm it actually still works post-move — it does (Steps 1–4 verified; Step 7's Claude
   API call was not exercised here since it requires a real `ANTHROPIC_API_KEY`, which this
   sandbox doesn't have).

3. **`patient_1.vcf` / `patient_1_report.json` are the synthetic demo patient**, not real
   patient data — confirmed by checking the VCF header and content directly (13 synthetic
   variants, no PHI). Safe to keep in `data/sample/` under version control, unlike anything
   that would later go in `data/real/` (gitignored by the project's own convention).

## Suggested next step

Decide whether `tnbc-genomics-agent`'s 10-kinase dictionary should be reconciled with the
90-kinase CTS panel (e.g., cross-check that the 10 genes here have consistent
inhibitor/resistance data against the DGIdb-sourced dictionary in `src/utils/`), or whether
it's deliberately staying a lighter-weight, independent tool. Either is defensible — just
worth being a decision rather than an unnoticed divergence.
