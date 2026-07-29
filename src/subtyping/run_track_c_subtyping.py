"""
run_track_c_subtyping.py

Runs the Track C TNBC subtyping pipeline against your REAL data:
    ~/rtk_nrtk_tnbc/data/processed/clinical.csv       (from assemble_tcga_brca_data.R)
    ~/rtk_nrtk_tnbc/data/processed/expression_tpm.csv (from assemble_tcga_brca_data.R)

Run from src/subtyping/ (needs src/data_loaders/gdc_local_data_loader.py and
src/regimen/cohort_wide_regimen_analysis.py's identify_tnbc_patients()
importable -- see the sys.path setup below):
    python3 run_track_c_subtyping.py
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).parent.parent / "data_loaders"))
sys.path.insert(0, str(Path(__file__).parent.parent / "regimen"))

from tnbc_de_novo_subtyping import run_tnbc_subtyping_pipeline, MARKER_GENE_SETS

try:
    from gdc_local_data_loader import load_r_exported_clinical
except ImportError:
    raise SystemExit(
        "Couldn't import gdc_local_data_loader.py -- confirm it's in src/data_loaders/."
    )

try:
    from cohort_wide_regimen_analysis import identify_tnbc_patients
except ImportError:
    raise SystemExit(
        "Couldn't import cohort_wide_regimen_analysis.py -- confirm it's in src/regimen/."
    )


def main():
    args = [a for a in sys.argv[1:] if a != "--full"]
    full_transcriptome = "--full" in sys.argv
    base_dir = args[0] if len(args) > 0 else str(Path.home() / "rtk_nrtk_tnbc")
    n_clusters = int(args[1]) if len(args) > 1 else 5

    # Real path confirmed on the actual machine: these live under a tcga_brca/ subfolder,
    # not directly in data/processed/ as originally guessed. A SEPARATE, differently-named
    # data/processed/tcga/clinical_data.csv also exists (the tcga vs tcga_brca ambiguity
    # flagged in assemble_tcga_brca_data.R's own comment and data/README.md) -- this uses
    # the tcga_brca/ one specifically, since its clinical.csv + expression_tpm.csv naming
    # matches assemble_tcga_brca_data.R's actual output convention exactly.
    clinical_path = str(Path(base_dir) / "data/processed/tcga_brca/clinical.csv")
    expression_filename = "expression_full_transcriptome.csv" if full_transcriptome else "expression_for_track_c.csv"
    expression_path = str(Path(base_dir) / "data/processed/tcga_brca" / expression_filename)

    print(f"Loading real clinical data from {clinical_path}...")
    clinical_df = load_r_exported_clinical(clinical_path)
    print(f"  {len(clinical_df)} patients in the full clinical file")

    # Real gap found and resolved on the actual machine: clinical.csv (from
    # GDCquery_clinic-style demographic/diagnosis/treatment data) has NO receptor-status
    # or PAM50 subtype columns at all -- confirmed by printing its real column list.
    # TCGA-BRCA molecular subtype/receptor status lives in a SEPARATE curated source
    # (TCGAbiolinks::TCGAquery_subtype(tumor="BRCA"), per the PanCanAtlas BRCA subtype
    # paper), with its own real column BRCA_Subtype_PAM50 -- use that file for TNBC
    # identification specifically, not clinical.csv.
    subtype_path = str(Path(base_dir) / "data/processed/tcga_brca/molecular_subtype.csv")

    print(f"\nIdentifying confirmed-TNBC patients (from {subtype_path}, BRCA_Subtype_PAM50)...")
    tnbc_patient_ids_raw = identify_tnbc_patients(subtype_path, subtype_column="BRCA_Subtype_PAM50")
    # Align to the same 12-character patient_id format load_r_exported_clinical()/
    # load_r_exported_expression() use, in case identify_tnbc_patients() returned full
    # sample-level barcodes rather than patient-level ones.
    tnbc_patient_ids = [pid[:12] for pid in tnbc_patient_ids_raw]
    print(f"  {len(tnbc_patient_ids)} confirmed-TNBC patients identified")

    n_overlap = len(set(tnbc_patient_ids) & set(clinical_df.index))
    if n_overlap < len(tnbc_patient_ids):
        print(f"  Note: only {n_overlap}/{len(tnbc_patient_ids)} TNBC patient IDs actually "
              f"found in clinical_df's index -- check for a patient-ID format mismatch "
              f"(e.g. barcode length, dashes vs. no dashes) if this gap is large.")

    print(f"\nLoading real expression data from {expression_path}...")
    # Real gap found and resolved: the R-based assembly script (assemble_tcga_brca_data.R)
    # crashed twice trying to build the full-transcriptome matrix via GDCprepare() (an OOM
    # issue confirmed to originate INSIDE GDCprepare() itself, not fixable by subsetting
    # genes afterward). Replaced with build_real_expression_for_track_c.py, which reads
    # the raw per-sample STAR-counts files directly and memory-safely (same fix your
    # kinase-scoring-pipeline README already documents for this exact GDCprepare() issue).
    # That script's output already has patient_id as the index (not a separate column),
    # so this loads it directly rather than via load_r_exported_expression()'s
    # R-specific column convention.
    expression_df = pd.read_csv(expression_path, index_col=0)
    print(f"  Expression matrix: {expression_df.shape[0]} patients x {expression_df.shape[1]} genes")

    print(f"\nRunning Track C subtyping pipeline (n_clusters={n_clusters}, "
          f"gene_selection_mode={'highly_variable' if full_transcriptome else 'marker_only'})...")
    result = run_tnbc_subtyping_pipeline(
        clinical_df, expression_df, tnbc_patient_ids, n_clusters=n_clusters,
        gene_selection_mode="highly_variable" if full_transcriptome else "marker_only",
    )

    print(f"\n=== Quality ===")
    print(f"Silhouette score: {result['quality']['silhouette_score']:.3f}")
    print(f"Cluster sizes: {result['quality']['cluster_sizes']}")

    print(f"\n=== Cluster characterization ===")
    print(result["cluster_characterization"].to_string())

    print(f"\n=== Patient assignments (first 20) ===")
    print(result["patient_assignments"].head(20).to_string())

    out_path = Path(base_dir) / "data/processed/tnbc_subtype_assignments.tsv"
    result["patient_assignments"].to_csv(out_path, sep="\t")
    print(f"\nWrote patient-level subtype assignments to {out_path}")

    print("\nHONEST REMINDER: no Lehmann/Burstein ground-truth label exists in this cohort to "
          "validate these calls against. Judge them by the silhouette score and marker "
          "enrichment table above, not by label alone. Try a few different n_clusters values "
          "(e.g. python3 run_track_c_subtyping.py ~/rtk_nrtk_tnbc 4  or  ...6) and compare "
          "silhouette scores before settling on one.")


if __name__ == "__main__":
    main()
