"""
build_real_expression_for_track_c.py

Replaces the R script's expression-assembly step entirely, after it crashed
twice (confirmed OOM-pattern: "100% then crash", consistent both times even
after subsetting genes post-GDCprepare(), meaning GDCprepare() itself -- not
anything downstream -- is what exhausts memory assembling the full ~1098 x
~60,000 SummarizedExperiment object).

This uses gdc_local_data_loader.py's build_expression_matrix(), which reads
each of your already-downloaded per-sample STAR-counts TSV files directly,
ONE AT A TIME, immediately restricted to just the ~126 needed genes (90
kinase panel + Track C's 36 markers) -- the same memory-safe approach your
tnbc-kinase-scoring-pipeline README already documents fixing this exact
GDCprepare() OOM crash for, just applied here to marker genes instead of
only the kinase panel. No R, no GDCprepare(), involved at all.

Run from src/data_loaders/:
    python3 build_real_expression_for_track_c.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from gdc_local_data_loader import build_expression_matrix


def main():
    args = [a for a in sys.argv[1:] if a != "--full"]
    base_dir = args[0] if len(args) > 0 else str(Path.home() / "rtk_nrtk_tnbc")
    full_transcriptome = "--full" in sys.argv

    root_path = str(Path(base_dir) / "data/raw/tcga_brca/GDCdata")
    sample_sheet_path = str(Path(base_dir) / "data/raw/tcga_brca/gdc_sample_sheet.tsv")

    if full_transcriptome:
        genes_needed = None
        print("Building FULL-TRANSCRIPTOME expression matrix (--full flag set) -- needed for "
              "gene_selection_mode='highly_variable' in Track C, since that mode needs a broad "
              "gene set to select the most-variable genes FROM. Confirmed memory-safe: this "
              "loader processes one file at a time (unlike GDCprepare(), which crashed twice) "
              "and a full ~1095 x ~20,000 matrix is only ~170MB as float64 -- nothing like the "
              "GDCprepare() OOM. Expect this to take noticeably longer than the restricted-gene "
              "run (more columns to parse per file), but not more memory-risky.\n")
    else:
        kinase_90_path = Path(base_dir) / "data/raw/kinases/kinase_90_list.txt"
        kinase_genes = []
        if kinase_90_path.exists():
            kinase_genes = [g.strip() for g in open(kinase_90_path) if g.strip()]
        else:
            print(f"Note: {kinase_90_path} not found -- proceeding with only Track C's marker genes.")

        track_c_marker_genes = [
            "AURKA", "AURKB", "PLK1", "CCNE1", "MYC", "CHEK1", "RAD51", "BRCA1",  # BL1
            "EGFR", "MET", "EPHA2", "PDGFRA", "KIT", "IGF1R",                      # BL2
            "VIM", "SNAI2", "TWIST1", "ZEB1", "CDH2", "FN1",                      # M
            "ALDH1A1", "ALDH1A3", "PDGFRB", "ABCB1",                              # MSL
            "AR", "FOXA1", "GATA3", "XBP1", "SPDEF",                              # LAR
            "CD3D", "CD8A", "CD274", "PDCD1", "CXCL9", "CXCL10", "GZMA",         # IM
        ]
        genes_needed = sorted(set(kinase_genes) | set(track_c_marker_genes))
        print(f"Building expression matrix restricted to {len(genes_needed)} genes "
              f"(90-kinase panel + Track C markers). Pass --full for the full transcriptome "
              f"instead (needed for gene_selection_mode='highly_variable').")

    print(f"Root path: {root_path}")
    print(f"Sample sheet: {sample_sheet_path}")
    print("Processing files one at a time (memory-safe) -- this may take a few minutes "
          "across ~1200+ files, but should never approach the memory ceiling that "
          "crashed GDCprepare().\n")

    matrix = build_expression_matrix(root_path, sample_sheet_path, gene_symbols=genes_needed)

    print(f"\nBuilt matrix: {matrix.shape[0]} patients x {matrix.shape[1]} genes")
    if genes_needed is not None:
        genes_actually_found = matrix.columns[matrix.notna().any()].tolist()
        genes_all_missing = matrix.columns[matrix.isna().all()].tolist()
        print(f"Genes with at least some real data: {len(genes_actually_found)}/{len(genes_needed)}")
        if genes_all_missing:
            print(f"Genes with NO data in any sample (symbol mismatch or genuinely absent): "
                  f"{genes_all_missing}")

    out_name = "expression_full_transcriptome.csv" if full_transcriptome else "expression_for_track_c.csv"
    out_path = Path(base_dir) / "data/processed/tcga_brca" / out_name
    matrix.to_csv(out_path)
    print(f"\nWrote real expression matrix to {out_path}")
    print("This is a SEPARATE file from expression_tpm.csv (the kinase-only file other "
          "parts of the pipeline already use) -- nothing existing was overwritten.")
    if not full_transcriptome:
        print("\nNext: point run_track_c_subtyping.py at this file instead of expression_tpm.csv.")
    else:
        print("\nNext: run_track_c_subtyping.py --full to use gene_selection_mode='highly_variable' "
              "against this file.")


if __name__ == "__main__":
    main()
