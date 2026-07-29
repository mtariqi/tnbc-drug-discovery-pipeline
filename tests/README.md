# tests/

Mirrors `src/` one subfolder at a time. Follows the convention already used throughout this
project (see the `_run_smoke_test()` functions inside `cohort_wide_regimen_analysis.py`,
`tnbc_specific_survival.py`, and `patient_level_survival_model.py`): synthetic data with a
**known** true structure, asserting the code recovers it — not just "does it run without
crashing."

Currently the smoke tests live inside each module's own `if __name__ == "__main__":` block
rather than as separate pytest files. That's fine for solo, real-data-verified development,
but before this repo has multiple contributors or CI, worth extracting into
`tests/<subfolder>/test_<module>.py` files pytest can discover and run as
`pytest tests/ -v` — a mechanical refactor, not a rewrite, since the actual test logic already
exists and already passes.
