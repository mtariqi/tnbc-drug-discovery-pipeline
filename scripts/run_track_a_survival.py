"""
Entry point: Track A patient-level survival model, full TCGA-BRCA cohort.

Wires src/data_loaders (once built) into src/survival/patient_level_survival_model.py.
Currently a template — the loader imports below will fail until
src/data_loaders/tcga_brca_loader.py exists (see its README for what to copy in).

Run from the repo root: python scripts/run_track_a_survival.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.survival.patient_level_survival_model import (
    build_patient_feature_matrix,
    train_and_evaluate_survival_model,
    stratify_by_risk,
    kaplan_meier_estimate,
    permutation_importance_cox,
)

# TODO once src/data_loaders/tcga_brca_loader.py exists:
# from src.data_loaders.tcga_brca_loader import load_clinical, load_expression
KINASE_PANEL_PATH = Path(__file__).resolve().parents[1] / "data" / "processed" / "kinase_panel.txt"


def main():
    if not KINASE_PANEL_PATH.exists():
        raise FileNotFoundError(
            f"{KINASE_PANEL_PATH} not found. Write your 90-kinase panel (one gene symbol per "
            f"line) there first, or point this script at wherever src/scoring/cts.py keeps it."
        )
    kinase_panel = KINASE_PANEL_PATH.read_text().split()

    # clinical_df = load_clinical(...)
    # expression_df = load_expression(...)
    # feature_df = build_patient_feature_matrix(clinical_df, expression_df, kinase_panel, "BRCA_Subtype_PAM50")
    # result = train_and_evaluate_survival_model(feature_df, clinical_df["duration"], clinical_df["event"])
    # print(f"Mean C-index: {result['mean_c_index']:.3f}")

    raise NotImplementedError(
        "Uncomment the lines above once src/data_loaders/tcga_brca_loader.py is built."
    )


if __name__ == "__main__":
    main()
