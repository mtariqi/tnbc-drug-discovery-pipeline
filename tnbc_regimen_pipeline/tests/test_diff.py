"""Tests for tnbc_regimen_pipeline.pipeline.diff."""

from __future__ import annotations

from tnbc_regimen_pipeline.pipeline.diff import diff_regimens_with_and_without_discovery


def test_diff_reports_hcos_delta_and_reason_for_change(
    real_maf_file, curated_gene_drugs, fake_resolve_tp53, fake_drug_graph_cls, fake_synergy_net_cls, fake_hcos, fake_mdcoe,
):
    merged_gene_drugs = {**curated_gene_drugs, "PTEN": ["alpelisib", "capivasertib"]}

    diff = diff_regimens_with_and_without_discovery(
        ["PATIENT-001", "PATIENT-002"], real_maf_file, ["EGFR", "PTEN", "TYK2"],
        curated_gene_drugs, merged_gene_drugs,
        fake_resolve_tp53, fake_drug_graph_cls, fake_synergy_net_cls, fake_hcos, fake_mdcoe,
    )

    assert list(diff["patient"]) == ["PATIENT-001"], "only the patient with PTEN altered should show a diff"
    row = diff.iloc[0]
    assert row["hcos_delta"] > 0, "HCOS should increase after adding a second PTEN drug"
    assert "capivasertib" in row["reason_for_change"]


def test_diff_returns_empty_frame_with_correct_columns_when_nothing_changes(
    real_maf_file, curated_gene_drugs, fake_resolve_tp53, fake_drug_graph_cls, fake_synergy_net_cls, fake_hcos, fake_mdcoe,
):
    diff = diff_regimens_with_and_without_discovery(
        ["PATIENT-001", "PATIENT-002"], real_maf_file, ["EGFR", "PTEN", "TYK2"],
        curated_gene_drugs, curated_gene_drugs,  # identical dict -> nothing should change
        fake_resolve_tp53, fake_drug_graph_cls, fake_synergy_net_cls, fake_hcos, fake_mdcoe,
    )
    assert diff.empty
    assert "hcos_delta" in diff.columns
    assert "reason_for_change" in diff.columns
