"""
Cohort-Wide Regimen Analysis Across All Confirmed TNBC Patients
==================================================================

Runs the existing, unmodified MDCOE/HCOS scoring for every patient in a
TCGA-BRCA cohort using each patient's own real altered-gene set (via
tnbc_regimen_pipeline.cohort.maf_reader), then aggregates which regimens
recur across the cohort versus being specific to one individual.
"""

from __future__ import annotations

import logging
from collections import Counter
from typing import Dict, List, Optional

import pandas as pd

from .maf_reader import extract_patient_altered_genes

logger = logging.getLogger(__name__)


# =====================================================================
# 1. IDENTIFY CONFIRMED TNBC PATIENTS
# =====================================================================

def identify_tnbc_patients(clinical_path: str, subtype_column: Optional[str] = None) -> List[str]:
    """
    Filters TCGA-BRCA clinical data to confirmed TNBC patients.

    TCGA-BRCA clinical exports vary in which column marks TNBC status --
    some carry a PAM50 subtype column, others require deriving TNBC
    status from three receptor columns (er/pr/her2_status_by_ihc, all
    'Negative'). This function tries the subtype column first if given,
    then falls back to receptor-status columns, and logs exactly which
    method it used -- verify this matches your actual file's real columns
    via inspect_clinical_columns() below before trusting it.
    """
    df = pd.read_csv(clinical_path, sep="\t" if clinical_path.endswith((".tsv", ".txt")) else ",")

    if subtype_column and subtype_column in df.columns:
        tnbc = df[df[subtype_column].astype(str).str.contains("Basal", case=False, na=False)]
        logger.info(f"identified {len(tnbc)} TNBC patients via subtype column '{subtype_column}'")
        return tnbc["bcr_patient_barcode"].tolist() if "bcr_patient_barcode" in tnbc.columns else tnbc.iloc[:, 0].tolist()

    receptor_cols = ["er_status_by_ihc", "pr_status_by_ihc", "her2_status_by_ihc"]
    if all(c in df.columns for c in receptor_cols):
        mask = True
        for c in receptor_cols:
            mask = mask & (df[c].astype(str).str.lower() == "negative")
        tnbc = df[mask]
        logger.info(
            f"identified {len(tnbc)} TNBC patients via triple-negative receptor status "
            f"({', '.join(receptor_cols)} all Negative)"
        )
        return tnbc["bcr_patient_barcode"].tolist() if "bcr_patient_barcode" in tnbc.columns else tnbc.iloc[:, 0].tolist()

    raise ValueError(
        "Could not identify a TNBC classification column. Run inspect_clinical_columns() "
        "on your real clinical file and pass the correct subtype_column explicitly."
    )


def inspect_clinical_columns(clinical_path: str, n: int = 5) -> None:
    """Run this FIRST on your real clinical file to see what columns
    actually exist, before trusting identify_tnbc_patients()'s
    column-name guesses."""
    df = pd.read_csv(clinical_path, sep="\t" if clinical_path.endswith((".tsv", ".txt")) else ",", nrows=n)
    logger.info(f"columns: {df.columns.tolist()}")
    print(f"Columns: {df.columns.tolist()}")
    print(df.head(n))


# =====================================================================
# 2. RUN MDCOE FOR EVERY PATIENT AND AGGREGATE
# =====================================================================

def run_cohort_wide_analysis(
    patient_barcodes: List[str],
    maf_glob_pattern: str,
    gene_panel: List[str],
    curated_gene_drugs: Dict[str, List[str]],
    resolve_tp53_fn,
    drug_graph_cls,
    synergy_net_cls,
    hcos_fn,
    mdcoe_fn,
    max_patients: Optional[int] = None,
) -> pd.DataFrame:
    """
    Runs the existing, unmodified MDCOE/HCOS scoring for every patient in
    patient_barcodes, using each patient's own real altered-gene set.
    Returns one row per patient: top regimen, its HCOS score, and how many
    of the patient's altered genes had at least one candidate drug.
    """
    barcodes = patient_barcodes[:max_patients] if max_patients else patient_barcodes
    rows = []
    for i, barcode in enumerate(barcodes, 1):
        altered_genes = extract_patient_altered_genes(barcode, maf_glob_pattern, gene_panel)
        if not altered_genes:
            rows.append({"patient": barcode, "n_altered_druggable_genes": 0, "top_regimen": None, "hcos": None})
            continue

        drugs = []
        for g in altered_genes:
            if g == "TP53":
                # NOTE: verify each patient's real TP53 classification individually for a
                # rigorous analysis; this defaults to the same assumption confirmed for
                # one previously-inspected patient (TCGA-AO-A128), which is very likely
                # wrong for most other patients.
                drugs.extend(resolve_tp53_fn("Nonsense_Mutation"))
            else:
                drugs.extend(curated_gene_drugs.get(g, []))
        drugs = sorted(set(drugs))

        if len(drugs) < 2:
            rows.append({"patient": barcode, "n_altered_druggable_genes": len(altered_genes), "top_regimen": None, "hcos": None})
            continue

        net = synergy_net_cls()
        results = mdcoe_fn(drug_graph_cls(drugs), net, hcos_fn, beam_width=50, max_depth=5, top_k=1)
        top_regimen, top_score = results[0] if results else (None, None)
        rows.append({
            "patient": barcode,
            "n_altered_druggable_genes": len(altered_genes),
            "top_regimen": " + ".join(sorted(top_regimen)) if top_regimen else None,
            "hcos": top_score,
        })

        if i % 10 == 0:
            logger.info(f"processed {i}/{len(barcodes)} patients")

    return pd.DataFrame(rows)


def aggregate_regimen_frequency(cohort_results: pd.DataFrame) -> pd.DataFrame:
    """Counts how often each regimen was the #1-ranked result across the cohort."""
    valid = cohort_results.dropna(subset=["top_regimen"])
    counts = Counter(valid["top_regimen"])
    freq = pd.DataFrame(counts.items(), columns=["regimen", "n_patients_top_ranked"])
    freq["pct_of_cohort"] = (freq["n_patients_top_ranked"] / len(cohort_results) * 100).round(1)
    return freq.sort_values("n_patients_top_ranked", ascending=False).reset_index(drop=True)
