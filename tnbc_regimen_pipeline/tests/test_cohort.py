"""Tests for tnbc_regimen_pipeline.cohort (MAF reading, cohort-wide
scoring, frequency aggregation)."""

from __future__ import annotations

from tnbc_regimen_pipeline.cohort import maf_reader, cohort_wide_regimen_analysis as cwra


def test_extract_patient_altered_genes_reads_real_gzipped_maf(real_maf_file):
    """Uses a REAL gzipped MAF file on disk (not mocked) to test the
    actual gzip + tab-delimited parsing logic."""
    genes_001 = maf_reader.extract_patient_altered_genes("PATIENT-001", real_maf_file, ["EGFR", "PTEN", "TYK2"])
    genes_002 = maf_reader.extract_patient_altered_genes("PATIENT-002", real_maf_file, ["EGFR", "PTEN", "TYK2"])
    assert genes_001 == ["EGFR", "PTEN"]
    assert genes_002 == ["EGFR"]


def test_extract_patient_altered_genes_respects_gene_panel(real_maf_file):
    """A gene altered for the patient but NOT in gene_panel must be excluded."""
    genes = maf_reader.extract_patient_altered_genes("PATIENT-001", real_maf_file, ["EGFR"])
    assert genes == ["EGFR"], "PTEN alteration exists for this patient but is outside the requested gene panel"


def test_run_cohort_wide_analysis_produces_expected_regimens(
    monkeypatch, real_maf_file, curated_gene_drugs, fake_resolve_tp53, fake_drug_graph_cls, fake_synergy_net_cls, fake_hcos, fake_mdcoe,
):
    result = cwra.run_cohort_wide_analysis(
        ["PATIENT-001", "PATIENT-002"], real_maf_file, ["EGFR", "PTEN", "TYK2"], curated_gene_drugs,
        fake_resolve_tp53, fake_drug_graph_cls, fake_synergy_net_cls, fake_hcos, fake_mdcoe,
    )
    row_001 = result[result["patient"] == "PATIENT-001"].iloc[0]
    row_002 = result[result["patient"] == "PATIENT-002"].iloc[0]

    assert row_001["n_altered_druggable_genes"] == 2
    assert set(row_001["top_regimen"].split(" + ")) == {"afatinib", "erlotinib", "alpelisib"}
    assert row_002["top_regimen"] == "afatinib + erlotinib"


def test_run_cohort_wide_analysis_handles_patient_with_no_druggable_alterations(
    tmp_path, curated_gene_drugs, fake_resolve_tp53, fake_drug_graph_cls, fake_synergy_net_cls, fake_hcos, fake_mdcoe,
):
    empty_maf = tmp_path / "empty.maf.gz"
    import gzip
    with gzip.open(empty_maf, "wt") as f:
        f.write("Hugo_Symbol\tVariant_Classification\tTumor_Sample_Barcode\n")

    result = cwra.run_cohort_wide_analysis(
        ["PATIENT-NONE"], str(empty_maf), ["EGFR", "PTEN"], curated_gene_drugs,
        fake_resolve_tp53, fake_drug_graph_cls, fake_synergy_net_cls, fake_hcos, fake_mdcoe,
    )
    row = result.iloc[0]
    assert row["n_altered_druggable_genes"] == 0
    assert row["top_regimen"] is None


def test_aggregate_regimen_frequency_counts_shared_regimens():
    import pandas as pd
    cohort_result = pd.DataFrame([
        {"patient": "P1", "n_altered_druggable_genes": 2, "top_regimen": "afatinib + alpelisib", "hcos": 0.2},
        {"patient": "P2", "n_altered_druggable_genes": 2, "top_regimen": "afatinib + alpelisib", "hcos": 0.2},
        {"patient": "P3", "n_altered_druggable_genes": 1, "top_regimen": None, "hcos": None},
    ])
    freq = cwra.aggregate_regimen_frequency(cohort_result)
    assert freq.iloc[0]["regimen"] == "afatinib + alpelisib"
    assert freq.iloc[0]["n_patients_top_ranked"] == 2
    assert freq.iloc[0]["pct_of_cohort"] == round(2 / 3 * 100, 1)
