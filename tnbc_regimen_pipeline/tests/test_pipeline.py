"""Tests for tnbc_regimen_pipeline.pipeline.full_pipeline (end-to-end
orchestration, exports, visual summary) and utils.validation."""

from __future__ import annotations

import json
import os
import xml.etree.ElementTree as ET

import pytest

from tnbc_regimen_pipeline.config import PipelineConfig
from tnbc_regimen_pipeline.pipeline.full_pipeline import run_full_pipeline
from tnbc_regimen_pipeline.pipeline.reproducibility import PIPELINE_VERSION
from tnbc_regimen_pipeline.utils.validation import validate_gene_drug_schema


def test_validate_gene_drug_schema_rejects_non_list_value():
    with pytest.raises(ValueError, match="must be a list"):
        validate_gene_drug_schema({"EGFR": "afatinib"})


def test_validate_gene_drug_schema_rejects_empty_drug_string():
    with pytest.raises(ValueError, match="non-string or empty"):
        validate_gene_drug_schema({"EGFR": ["afatinib", ""]})


def test_validate_gene_drug_schema_accepts_well_formed_dict():
    validate_gene_drug_schema({"EGFR": ["afatinib"], "PTEN": []})  # should not raise


def test_run_full_pipeline_end_to_end(
    tmp_path, monkeypatch, real_maf_file, curated_gene_drugs, real_gene_drugs,
    fake_resolve_tp53, fake_drug_graph_cls, fake_synergy_net_cls, fake_hcos, fake_mdcoe,
    mock_requests_get, mock_requests_post,
):
    from tnbc_regimen_pipeline.discovery import agentic_regimen_discovery as ard
    monkeypatch.setattr(ard.requests, "get", mock_requests_get)
    monkeypatch.setattr(ard.requests, "post", mock_requests_post)

    config = PipelineConfig(max_papers_per_gene=5, output_dir=str(tmp_path / "outputs"), max_workers=4)

    result = run_full_pipeline(
        genes=["EGFR", "PTEN", "TYK2"],
        curated_gene_drugs=curated_gene_drugs,
        real_gene_drugs=real_gene_drugs,
        patient_barcodes=["PATIENT-001", "PATIENT-002"],
        maf_glob_pattern=real_maf_file,
        gene_panel=["EGFR", "PTEN", "TYK2"],
        resolve_tp53_fn=fake_resolve_tp53,
        drug_graph_cls=fake_drug_graph_cls,
        synergy_net_cls=fake_synergy_net_cls,
        hcos_fn=fake_hcos,
        mdcoe_fn=fake_mdcoe,
        config=config,
    )

    # discovery correctly isolated the real gap gene
    assert set(result["discovery_result"]["gene"].unique()) == {"PTEN"}

    # provenance correctly attributes the PMID; curated drugs stay sourceless
    assert result["provenance"]["PTEN"]["capivasertib"] == ["999"]
    assert result["provenance"]["PTEN"]["alpelisib"] == []

    # the patient with PTEN altered picked up the discovered, confirmed drug
    p1_regimen = result["cohort_result"].loc[result["cohort_result"]["patient"] == "PATIENT-001", "top_regimen"].iloc[0]
    assert "capivasertib" in p1_regimen

    # diff isolates exactly that patient
    assert list(result["diff"]["patient"]) == ["PATIENT-001"]

    # exports landed on disk
    for name in ["discovery_result", "cohort_result", "regimen_frequency", "regimen_diff", "merged_gene_drugs", "provenance"]:
        assert name in result["export_paths"]
        assert os.path.exists(result["export_paths"][name])

    # exported JSON carries version-stamping metadata
    with open(result["export_paths"]["provenance"]) as f:
        exported = json.load(f)
    assert exported["_metadata"]["pipeline_version"] == PIPELINE_VERSION
    assert "git_commit" in exported["_metadata"]
    assert "run_timestamp_utc" in exported["_metadata"]

    # visuals exist and the SVG is well-formed XML
    assert os.path.exists(result["visual_paths"]["bar_chart"])
    assert os.path.exists(result["visual_paths"]["flow_svg"])
    ET.parse(result["visual_paths"]["flow_svg"])  # raises if malformed


def test_run_full_pipeline_rejects_malformed_curated_dict(
    tmp_path, real_maf_file, real_gene_drugs,
    fake_resolve_tp53, fake_drug_graph_cls, fake_synergy_net_cls, fake_hcos, fake_mdcoe,
):
    config = PipelineConfig(output_dir=str(tmp_path / "outputs"))
    with pytest.raises(ValueError):
        run_full_pipeline(
            genes=["EGFR"],
            curated_gene_drugs={"EGFR": "afatinib"},  # malformed: string, not list
            real_gene_drugs=real_gene_drugs,
            patient_barcodes=["PATIENT-001"],
            maf_glob_pattern=real_maf_file,
            gene_panel=["EGFR"],
            resolve_tp53_fn=fake_resolve_tp53,
            drug_graph_cls=fake_drug_graph_cls,
            synergy_net_cls=fake_synergy_net_cls,
            hcos_fn=fake_hcos,
            mdcoe_fn=fake_mdcoe,
            config=config,
        )
