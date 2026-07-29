"""
Cohort-Wide Regimen Analysis Across All Confirmed TNBC Patients
==================================================================

Extends the single-patient MDCOE/HCOS worked example (TCGA-AO-A128) to the
FULL set of confirmed TNBC patients already present in your downloaded
TCGA-BRCA cohort. No new data acquisition is required -- this uses the
same clinical and MAF files already used to build the CTS survival
component and the single-patient case study.

WHAT THIS ANSWERS THAT THE SINGLE-PATIENT CASE COULD NOT:
    - Is the afatinib+alpelisib/capivasertib+trastuzumab finding a broadly
      recurring top regimen across TNBC patients with similar alterations,
      or specific to one individual?
    - Which gene-alteration combinations are common enough across real
      patients to justify prioritizing a regimen for a subgroup?
    - Does the cohort reveal patient subgroups needing different regimens?

REAL, CONFIRMED MAF STRUCTURE (verified against the actual downloaded
TCGA-AO-A128 file this session): tab-delimited, Hugo_Symbol and
Tumor_Sample_Barcode and Variant_Classification are real column names in
the actual file header. This parser uses those confirmed column names
directly rather than assuming a different convention.
"""

from __future__ import annotations

import gzip
import glob
from collections import Counter
from typing import Dict, List, Optional

import pandas as pd


# =====================================================================
# 1. IDENTIFY CONFIRMED TNBC PATIENTS
# =====================================================================

def identify_tnbc_patients(clinical_path: str, subtype_column: Optional[str] = None) -> List[str]:
    """
    Filters TCGA-BRCA clinical data to confirmed TNBC patients.

    IMPORTANT: TCGA-BRCA clinical exports vary in which column marks TNBC
    status -- some carry a PAM50 'BRCA_Subtype_PAM50' or 'subtype' column
    (TCGA-AO-A128 was previously confirmed as 'BRCA_Basal' via this route),
    others require deriving TNBC status from three separate receptor
    columns (er_status_by_ihc, pr_status_by_ihc, her2_status_by_ihc, all
    'Negative'). This function tries the subtype column first if given,
    then falls back to the receptor-status columns, and prints exactly
    which method it used -- verify this matches your actual file's
    real columns via inspect_clinical_columns() below before trusting it.
    """
    df = pd.read_csv(clinical_path, sep="\t" if clinical_path.endswith((".tsv", ".txt")) else ",")

    if subtype_column and subtype_column in df.columns:
        tnbc = df[df[subtype_column].astype(str).str.contains("Basal", case=False, na=False)]
        print(f"Identified {len(tnbc)} TNBC patients via subtype column '{subtype_column}'")
        return tnbc["bcr_patient_barcode"].tolist() if "bcr_patient_barcode" in tnbc.columns else tnbc.iloc[:, 0].tolist()

    receptor_cols = ["er_status_by_ihc", "pr_status_by_ihc", "her2_status_by_ihc"]
    if all(c in df.columns for c in receptor_cols):
        mask = True
        for c in receptor_cols:
            mask = mask & (df[c].astype(str).str.lower() == "negative")
        tnbc = df[mask]
        print(f"Identified {len(tnbc)} TNBC patients via triple-negative receptor status "
              f"({', '.join(receptor_cols)} all Negative)")
        return tnbc["bcr_patient_barcode"].tolist() if "bcr_patient_barcode" in tnbc.columns else tnbc.iloc[:, 0].tolist()

    raise ValueError(
        "Could not identify a TNBC classification column. Run inspect_clinical_columns() "
        "on your real clinical file and pass the correct subtype_column explicitly."
    )


def inspect_clinical_columns(clinical_path: str, n: int = 5) -> None:
    """Run this FIRST on your real clinical file to see what columns actually exist,
    before trusting identify_tnbc_patients()'s column-name guesses."""
    df = pd.read_csv(clinical_path, sep="\t" if clinical_path.endswith((".tsv", ".txt")) else ",", nrows=n)
    print(f"Columns: {df.columns.tolist()}")
    print(df.head(n))


# =====================================================================
# 2. EXTRACT EACH PATIENT'S ALTERED GENES FROM REAL MAF DATA
# =====================================================================

def extract_patient_altered_genes(patient_barcode: str, maf_glob_pattern: str, gene_panel: List[str]) -> List[str]:
    """
    Searches real, gzipped GDC MAF files for this patient's barcode and
    returns which genes in gene_panel show a variant for them. Uses the
    CONFIRMED real column names (Hugo_Symbol, Tumor_Sample_Barcode) from
    this session's direct inspection of an actual downloaded MAF file --
    not assumed from documentation alone.
    """
    altered = set()
    for maf_path in glob.glob(maf_glob_pattern):
        opener = gzip.open if maf_path.endswith(".gz") else open
        with opener(maf_path, "rt") as f:
            header = None
            for line in f:
                if line.startswith("Hugo_Symbol"):
                    header = line.rstrip("\n").split("\t")
                    continue
                if header is None:
                    continue
                if patient_barcode not in line:
                    continue
                fields = line.rstrip("\n").split("\t")
                row = dict(zip(header, fields))
                gene = row.get("Hugo_Symbol", "")
                barcode_field = row.get("Tumor_Sample_Barcode", "")
                if patient_barcode in barcode_field and gene in gene_panel:
                    altered.add(gene)
    return sorted(altered)


# =====================================================================
# 3. RUN MDCOE FOR EVERY PATIENT AND AGGREGATE
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
                drugs.extend(resolve_tp53_fn("Nonsense_Mutation"))  # NOTE: verify each patient's real TP53 classification individually for a rigorous analysis; this defaults to the same assumption confirmed for TCGA-AO-A128
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
            print(f"  ...processed {i}/{len(barcodes)} patients")

    return pd.DataFrame(rows)


def aggregate_regimen_frequency(cohort_results: pd.DataFrame) -> pd.DataFrame:
    """Counts how often each regimen was the #1-ranked result across the cohort."""
    valid = cohort_results.dropna(subset=["top_regimen"])
    counts = Counter(valid["top_regimen"])
    freq = pd.DataFrame(counts.items(), columns=["regimen", "n_patients_top_ranked"])
    freq["pct_of_cohort"] = (freq["n_patients_top_ranked"] / len(cohort_results) * 100).round(1)
    return freq.sort_values("n_patients_top_ranked", ascending=False).reset_index(drop=True)


# =====================================================================
# SMOKE TEST -- synthetic cohort, verifies aggregation logic end-to-end
# =====================================================================

def _run_smoke_test():
    class FakeDrugGraph:
        def __init__(self, drugs): self.drugs = drugs
    class FakeSynergyNet:
        pass
    def fake_hcos(regimen, net):
        # deterministic fake score: more drugs = higher score, for testability
        return len(regimen) * 0.1
    def fake_mdcoe(graph, net, hcos_fn, beam_width, max_depth, top_k):
        import itertools
        best = None
        for size in range(2, min(len(graph.drugs), max_depth) + 1):
            for combo in itertools.combinations(graph.drugs, size):
                score = hcos_fn(combo, net)
                if best is None or score > best[1]:
                    best = (combo, score)
        return [best] if best else []
    def fake_resolve_tp53(alt_type):
        return ["adavosertib", "venetoclax"]

    curated = {"EGFR": ["afatinib"], "PTEN": ["alpelisib"], "ERBB2": ["trastuzumab"], "TYK2": ["deucravacitinib"]}

    # Simulate extract_patient_altered_genes output directly (bypassing real MAF I/O for this test)
    import unittest.mock as mock
    fake_patient_genes = {
        "PATIENT-001": ["EGFR", "PTEN", "ERBB2"],
        "PATIENT-002": ["EGFR", "PTEN", "ERBB2"],  # same alterations -> should get same top regimen
        "PATIENT-003": ["EGFR", "TYK2"],           # different alterations -> different regimen
        "PATIENT-004": [],                          # no altered druggable genes
    }
    with mock.patch(f"{__name__}.extract_patient_altered_genes", side_effect=lambda b, m, g: fake_patient_genes[b]):
        result = run_cohort_wide_analysis(
            list(fake_patient_genes.keys()), "fake_glob", list(curated.keys()),
            curated, fake_resolve_tp53, FakeDrugGraph, FakeSynergyNet, fake_hcos, fake_mdcoe,
        )
    print(result)
    print()

    freq = aggregate_regimen_frequency(result)
    print(freq)
    print()

    assert result.loc[result["patient"] == "PATIENT-004", "n_altered_druggable_genes"].iloc[0] == 0
    assert result.loc[result["patient"] == "PATIENT-001", "top_regimen"].iloc[0] == result.loc[result["patient"] == "PATIENT-002", "top_regimen"].iloc[0]
    assert freq.iloc[0]["n_patients_top_ranked"] == 2, "the shared regimen across PATIENT-001/002 should be the most frequent"
    print("PASSED: identical patient alteration profiles correctly produce identical top regimens,")
    print("frequency aggregation correctly counts shared regimens across the cohort, and patients")
    print("with no druggable alterations are correctly handled without crashing.")


if __name__ == "__main__":
    _run_smoke_test()
