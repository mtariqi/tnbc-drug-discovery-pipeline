"""
TNBC-Specific Survival Refinement for CTS
============================================

The original CTS survival component pools ALL TCGA-BRCA patients
(~1095-1098, spanning every breast cancer subtype) to compute each
kinase's Cox hazard ratio and log-rank association with survival. This
is a real, defensible refinement: recompute that same association using
ONLY the confirmed-TNBC subset (~150-200 real patients, restricted via
identify_tnbc_patients() in cohort_wide_regimen_analysis.py), and compare
against the original pooled-cohort result.

WHY THIS MATTERS, HONESTLY:
    A kinase's survival association in the full, mixed BRCA cohort could
    be driven mostly by hormone-receptor-positive or HER2+ patients, who
    make up the large majority of that cohort -- diluting or masking a
    real, different association specific to TNBC. Restricting to the
    confirmed-TNBC subset directly tests whether this happens, using real
    patient-level expression and survival data already downloaded for
    this project. This does NOT invent a drug-response label that
    doesn't exist in TCGA (a genuine data limitation, stated plainly) --
    it uses the real survival/expression signal TCGA actually contains,
    restricted to the population this whole project is about.

This module deliberately reuses run_survival_pipeline() from
tcga_brca_survival_pipeline.py UNMODIFIED -- the Cox/log-rank/RMST
machinery itself is already validated; what's new here is only the
TNBC-restriction and full-vs-TNBC comparison logic.
"""

from __future__ import annotations

from typing import Dict, List

import pandas as pd


# =====================================================================
# 1. RESTRICT COHORT TO CONFIRMED TNBC PATIENTS
# =====================================================================

def restrict_to_tnbc(
    clinical_df: pd.DataFrame,
    expression_df: pd.DataFrame,
    tnbc_patient_ids: List[str],
) -> tuple:
    """
    Filters both the clinical and expression DataFrames down to only the
    confirmed-TNBC patients (from identify_tnbc_patients()). Both inputs
    are assumed indexed or joinable on patient_id, matching
    run_survival_pipeline()'s existing expectations.
    """
    clinical_restricted = clinical_df[clinical_df.index.isin(tnbc_patient_ids)].copy()
    expression_restricted = expression_df[expression_df.index.isin(tnbc_patient_ids)].copy()

    n_found = len(clinical_restricted)
    n_requested = len(tnbc_patient_ids)
    if n_found < n_requested:
        print(f"Note: {n_found}/{n_requested} requested TNBC patient IDs found in the clinical data "
              f"(the rest may be missing survival/expression records, or an ID-format mismatch -- "
              f"worth checking if this gap is large).")

    return clinical_restricted, expression_restricted


# =====================================================================
# 2. RUN THE EXISTING, VALIDATED SURVIVAL PIPELINE ON THE TNBC SUBSET
# =====================================================================

def run_tnbc_specific_survival_analysis(
    clinical_df: pd.DataFrame,
    expression_df: pd.DataFrame,
    genes: List[str],
    tnbc_patient_ids: List[str],
    run_survival_pipeline_fn,
    adjust_for=("stage_numeric", "age"),
) -> pd.DataFrame:
    """
    Restricts to the confirmed-TNBC subset, then calls the existing,
    unmodified run_survival_pipeline() (pass the real function imported
    from tcga_brca_survival_pipeline.py). Returns the same column
    structure as the original, full-cohort CTS survival computation, so
    it can be directly compared or substituted in.
    """
    clinical_tnbc, expression_tnbc = restrict_to_tnbc(clinical_df, expression_df, tnbc_patient_ids)
    print(f"Running survival analysis on {len(clinical_tnbc)} confirmed-TNBC patients "
          f"(vs. the full ~1095-1098-patient pooled BRCA cohort used in the original CTS computation)")
    return run_survival_pipeline_fn(clinical_tnbc, expression_tnbc, genes, adjust_for=adjust_for)


# =====================================================================
# 3. COMPARE FULL-COHORT VS. TNBC-RESTRICTED RESULTS
# =====================================================================

def compare_full_vs_tnbc_survival(
    full_cohort_result: pd.DataFrame,
    tnbc_result: pd.DataFrame,
    hr_direction_flip_threshold: float = 1.0,
) -> pd.DataFrame:
    """
    Compares each kinase's Cox hazard ratio between the full pooled-BRCA
    result (already computed for the original CTS) and the new
    TNBC-restricted result. Flags kinases where the association's
    DIRECTION flips (HR crosses 1.0 -- protective in one population,
    harmful in the other) as the most scientifically interesting cases,
    since a magnitude-only difference is expected from sample-size noise
    alone, while a direction flip is a stronger, more specific signal.
    """
    joined = full_cohort_result[["cox_hr"]].rename(columns={"cox_hr": "cox_hr_full_cohort"}).join(
        tnbc_result[["cox_hr"]].rename(columns={"cox_hr": "cox_hr_tnbc_only"}), how="inner"
    )
    joined["direction_flipped"] = (
        (joined["cox_hr_full_cohort"] - hr_direction_flip_threshold)
        * (joined["cox_hr_tnbc_only"] - hr_direction_flip_threshold)
    ) < 0
    joined["abs_hr_difference"] = (joined["cox_hr_full_cohort"] - joined["cox_hr_tnbc_only"]).abs()
    return joined.sort_values("direction_flipped", ascending=False)


# =====================================================================
# SMOKE TEST -- verifies filtering + comparison logic without needing
# a real lifelines install (mocks the already-validated pipeline call)
# =====================================================================

def _run_smoke_test():
    # Synthetic full cohort: 10 patients, 3 TNBC + 7 non-TNBC
    clinical_df = pd.DataFrame(index=[f"P{i:03d}" for i in range(10)])
    expression_df = pd.DataFrame({"EGFR": range(10)}, index=[f"P{i:03d}" for i in range(10)])
    tnbc_ids = ["P000", "P001", "P002"]

    print("=== Testing restrict_to_tnbc() ===")
    clin_restricted, expr_restricted = restrict_to_tnbc(clinical_df, expression_df, tnbc_ids)
    print(f"Restricted clinical: {len(clin_restricted)} patients, restricted expression: {len(expr_restricted)} patients")
    assert len(clin_restricted) == 3
    assert set(clin_restricted.index) == set(tnbc_ids)
    print("PASSED: correctly restricted to only the 3 TNBC patients.\n")

    print("=== Testing restrict_to_tnbc() with a missing-patient gap ===")
    tnbc_ids_with_gap = tnbc_ids + ["P999_NOT_IN_COHORT"]
    clin_restricted2, _ = restrict_to_tnbc(clinical_df, expression_df, tnbc_ids_with_gap)
    assert len(clin_restricted2) == 3, "should still correctly find the 3 real patients despite 1 missing ID"
    print("PASSED: correctly handled a requested patient ID not present in the data, without crashing.\n")

    print("=== Testing run_tnbc_specific_survival_analysis() with a mocked pipeline call ===")
    def fake_run_survival_pipeline(clinical, expression, genes, adjust_for):
        # Confirms the function actually received the RESTRICTED (3-patient) data, not the full 10
        assert len(clinical) == 3, f"expected the restricted 3-patient cohort, got {len(clinical)}"
        return pd.DataFrame({"cox_hr": [1.5]}, index=["EGFR"])

    result = run_tnbc_specific_survival_analysis(
        clinical_df, expression_df, ["EGFR"], tnbc_ids, fake_run_survival_pipeline,
    )
    print(result)
    print("PASSED: correctly passed the TNBC-restricted subset (not the full cohort) into the survival pipeline.\n")

    print("=== Testing compare_full_vs_tnbc_survival() direction-flip detection ===")
    full_result = pd.DataFrame({"cox_hr": [1.8, 0.6, 1.1]}, index=["EGFR", "PTK2", "STABLE_GENE"])
    tnbc_result_df = pd.DataFrame({"cox_hr": [0.5, 0.55, 1.15]}, index=["EGFR", "PTK2", "STABLE_GENE"])
    # EGFR: 1.8 -> 0.5, crosses 1.0 -- direction flip (harmful in full cohort, protective in TNBC)
    # PTK2: 0.6 -> 0.55, stays below 1.0 -- no flip, just magnitude change
    # STABLE_GENE: 1.1 -> 1.15, stays above 1.0 -- no flip
    comparison = compare_full_vs_tnbc_survival(full_result, tnbc_result_df)
    print(comparison)
    assert comparison.loc["EGFR", "direction_flipped"] == True
    assert comparison.loc["PTK2", "direction_flipped"] == False
    assert comparison.loc["STABLE_GENE", "direction_flipped"] == False
    print("\nPASSED: correctly identifies EGFR's direction flip (HR crosses 1.0) as the scientifically")
    print("interesting case, while PTK2's magnitude-only shift is correctly NOT flagged as a flip.")


if __name__ == "__main__":
    _run_smoke_test()
