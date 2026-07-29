# =====================================================================
# Assemble TCGA-BRCA expression + clinical + mutation data into CSVs
# for the Python kinase/survival pipeline. Reuses your already-downloaded
# 21.4GB (no re-download) -- just re-fetches the lightweight query
# metadata and clinical/mutation data.
# =====================================================================

library(TCGAbiolinks)
library(SummarizedExperiment)

# !! CONFIRM this matches your real download path (see chat -- your
# screenshot showed .../raw/tcga_brca/... but your script used .../raw/tcga/...)
data_dir <- "~/rtk_nrtk_tnbc/data/raw/tcga_brca/GDCdata"
out_dir <- "~/rtk_nrtk_tnbc/data/processed"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------
# 1. EXPRESSION -- rebuild the exact same query (does NOT re-download;
#    GDCquery() only fetches lightweight metadata), then GDCprepare()
#    reads your already-downloaded files.
# ---------------------------------------------------------------------
expr_query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

expr_data <- tryCatch(
  {
    GDCprepare(expr_query, directory = data_dir, summarizedExperiment = TRUE)
  },
  error = function(e) {
    cat("GDCprepare() failed with:\n", conditionMessage(e), "\n")
    cat("This is a known issue on some TCGAbiolinks/readr version combos.\n")
    cat("Try: install.packages('readr') to update, or:\n")
    cat("BiocManager::install('TCGAbiolinks', update = TRUE, ask = FALSE)\n")
    stop("Fix the above, then re-run this script.")
  }
)

cat("Available assays:", paste(assayNames(expr_data), collapse = ", "), "\n")

# tpm_unstrand (NOT tpm_unstranded -- that's the raw per-file column name;
# GDCprepare renames it without the trailing '-ed' in the assembled object)
tpm_matrix <- assay(expr_data, "tpm_unstrand")

# Row names are Ensembl IDs with version suffix (e.g. 'ENSG00000146648.20').
# Attach gene symbols from rowData so the Python side can index by symbol.
gene_symbols <- rowData(expr_data)$gene_name
rownames(tpm_matrix) <- make.unique(gene_symbols)

# ---------------------------------------------------------------------
# REAL FIX for the OOM crash: subset to only the genes actually needed
# (90-kinase panel + Track C's 36 marker genes) BEFORE the expensive
# transpose + as.data.frame() conversion below. The original script
# materialized the FULL ~20,000-gene x ~1098-sample matrix as a data.frame,
# which is almost certainly what crashed the R session (confirmed pattern:
# progress reached 100%, i.e. GDCprepare() succeeded, then died -- exactly
# where this data.frame conversion happens). Subsetting first means this
# step only ever has to handle a few hundred genes, not tens of thousands.
#
# NOTE: I could not run/test this R script myself (no R available in the
# environment I'm working in) -- unlike every Python change in this
# project, this has NOT been execute-verified. Check the printed gene-count
# and "genes not found" messages below carefully before trusting the output.
# ---------------------------------------------------------------------
kinase_90_path <- "~/rtk_nrtk_tnbc/data/raw/kinases/kinase_90_list.txt"
kinase_genes <- if (file.exists(path.expand(kinase_90_path))) {
  readLines(path.expand(kinase_90_path))
} else {
  cat("Warning:", kinase_90_path, "not found -- proceeding with only the Track C marker genes.\n")
  character(0)
}

track_c_marker_genes <- c(
  "AURKA", "AURKB", "PLK1", "CCNE1", "MYC", "CHEK1", "RAD51", "BRCA1",  # BL1
  "EGFR", "MET", "EPHA2", "PDGFRA", "KIT", "IGF1R",                     # BL2
  "VIM", "SNAI2", "TWIST1", "ZEB1", "CDH2", "FN1",                      # M
  "ALDH1A1", "ALDH1A3", "PDGFRB", "ABCB1",                              # MSL (VIM/IGF1R already listed above)
  "AR", "FOXA1", "GATA3", "XBP1", "SPDEF",                              # LAR
  "CD3D", "CD8A", "CD274", "PDCD1", "CXCL9", "CXCL10", "GZMA"           # IM
)

genes_needed <- unique(c(kinase_genes, track_c_marker_genes))
genes_found <- intersect(genes_needed, rownames(tpm_matrix))
genes_missing <- setdiff(genes_needed, rownames(tpm_matrix))

cat("Subsetting expression matrix to", length(genes_found), "/", length(genes_needed),
    "needed genes (90-kinase panel + Track C markers), to avoid materializing the full\n")
cat("~", nrow(tpm_matrix), "-gene matrix as a data.frame (the likely cause of the earlier crash).\n")
if (length(genes_missing) > 0) {
  cat("Note:", length(genes_missing), "genes not found in this expression data (symbol mismatch",
      "or genuinely absent):", paste(genes_missing, collapse = ", "), "\n")
}

tpm_matrix <- tpm_matrix[genes_found, , drop = FALSE]

# Transpose so patients are rows (matches build_expression_matrix()'s output shape)
tpm_df <- as.data.frame(t(tpm_matrix))
tpm_df$patient_id <- substr(rownames(tpm_df), 1, 12)  # TCGA barcode -> patient id
write.csv(tpm_df, file.path(out_dir, "expression_tpm.csv"), row.names = TRUE)
cat("Wrote", nrow(tpm_df), "samples x", ncol(tpm_df) - 1, "genes to expression_tpm.csv\n")

# ---------------------------------------------------------------------
# 2. CLINICAL -- dedicated query (small download, not part of the 21.4GB;
#    this is the authoritative source, more complete than GDCprepare's
#    auto-attached colData).
# ---------------------------------------------------------------------
clinical_df <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")

cat("\nClinical columns available (check these match what")
cat(" derive_survival_time_event() expects):\n")
print(intersect(
  colnames(clinical_df),
  c("vital_status", "days_to_death", "days_to_last_follow_up",
    "ajcc_pathologic_stage", "age_at_index", "age_at_diagnosis")
))

write.csv(clinical_df, file.path(out_dir, "clinical.csv"), row.names = FALSE)
cat("Wrote", nrow(clinical_df), "patients to clinical.csv\n")

# ---------------------------------------------------------------------
# 3. MUTATIONS -- separate query + download (NOT part of your original
#    21.4GB expression-only download). MAF files are much smaller than
#    RNA-seq data, this should be quick.
# ---------------------------------------------------------------------
maf_query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  access = "open"
)
GDCdownload(maf_query, directory = data_dir)
maf_data <- GDCprepare(maf_query, directory = data_dir)

write.csv(maf_data, file.path(out_dir, "mutations_maf.csv"), row.names = FALSE)
cat("Wrote", nrow(maf_data), "mutation records to mutations_maf.csv\n")

cat("\nDone. Three files in", out_dir, ":\n")
cat("  expression_tpm.csv  -> patient x gene TPM matrix\n")
cat("  clinical.csv        -> vital_status, days_to_death, stage, age, etc.\n")
cat("  mutations_maf.csv   -> long mutation table (Hugo_Symbol, Tumor_Sample_Barcode, ...)\n")
cat("These three CSVs are what the Python pipeline (tcga_brca_survival_pipeline.py,\n")
cat("kinase_scoring_pipeline.py) expects -- load them with pandas.read_csv().\n")
