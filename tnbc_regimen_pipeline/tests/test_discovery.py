"""Tests for tnbc_regimen_pipeline.discovery (extraction, search, DGIdb
confirmation, caching, parallelization, provenance)."""

from __future__ import annotations

import shutil
import tempfile

import pandas as pd
import pytest

from tnbc_regimen_pipeline.discovery import agentic_regimen_discovery as ard
from tnbc_regimen_pipeline.discovery import caching, parallelization, provenance


def test_extract_candidate_drugs_from_abstract_excludes_generic_terms():
    abstract = (
        "In this study, we evaluated the combination of alpelisib and "
        "trastuzumab in PIK3CA-mutant, HER2-low breast cancer models. "
        "Patients were also treated with standard chemotherapy and pain "
        "medication."
    )
    candidates = ard.extract_candidate_drugs_from_abstract(abstract)
    assert "alpelisib" in candidates
    assert "trastuzumab" in candidates
    assert "chemotherapy" not in candidates, "generic terms must never be extracted as drug candidates"


def test_identify_coverage_gaps(curated_gene_drugs, real_gene_drugs):
    gaps = ard.identify_coverage_gaps(
        ["EGFR", "PTEN", "TYK2"], curated_gene_drugs, real_gene_drugs, min_candidates=2
    )
    assert gaps == ["PTEN"], "only PTEN has fewer than 2 total candidates across curated+real"


def test_search_pubmed_for_combination_therapy(monkeypatch, mock_requests_get):
    monkeypatch.setattr(ard.requests, "get", mock_requests_get)
    papers = ard.search_pubmed_for_combination_therapy("PTEN", max_results=2)
    assert len(papers) == 1
    assert papers[0]["pmid"] == "999"
    assert "alpelisib" in papers[0]["abstract"].lower()


def test_confirm_candidate_via_dgidb(monkeypatch, mock_requests_post):
    monkeypatch.setattr(ard.requests, "post", mock_requests_post)
    assert ard.confirm_candidate_via_dgidb("capivasertib", "PTEN") is True
    assert ard.confirm_candidate_via_dgidb("totally_made_up_drug", "PTEN") is False


def test_merge_into_candidate_pool_only_merges_confirmed():
    discovery_result = pd.DataFrame([
        {"gene": "PTEN", "candidate_drug": "capivasertib", "pmid": "111", "dgidb_confirmed": True, "abstract_snippet": "..."},
        {"gene": "PTEN", "candidate_drug": "fake_unconfirmed_drug", "pmid": "222", "dgidb_confirmed": False, "abstract_snippet": "..."},
    ])
    merged = ard.merge_into_candidate_pool(discovery_result, {"PTEN": ["alpelisib"]})
    assert "capivasertib" in merged["PTEN"]
    assert "fake_unconfirmed_drug" not in merged["PTEN"]


def test_caching_avoids_duplicate_network_calls(monkeypatch, mock_requests_get):
    call_count = {"n": 0}

    def counting_get(url, params=None, timeout=None):
        call_count["n"] += 1
        return mock_requests_get(url, params=params, timeout=timeout)

    monkeypatch.setattr(ard.requests, "get", counting_get)

    tmp_dir = tempfile.mkdtemp()
    try:
        caching.install_cache(cache_dir=tmp_dir)
        first = ard.search_pubmed_for_combination_therapy("PTEN", max_results=2)
        second = ard.search_pubmed_for_combination_therapy("PTEN", max_results=2)
        assert first == second
        assert call_count["n"] == 2, "second identical call should be served from cache (esearch+efetch counted once)"
    finally:
        caching.uninstall_cache()
        shutil.rmtree(tmp_dir, ignore_errors=True)


def test_parallel_discovery_only_searches_gap_genes(monkeypatch, mock_requests_get, mock_requests_post, curated_gene_drugs, real_gene_drugs):
    monkeypatch.setattr(ard.requests, "get", mock_requests_get)
    monkeypatch.setattr(ard.requests, "post", mock_requests_post)

    result = parallelization.parallel_discovery(
        ["EGFR", "PTEN", "TYK2"], curated_gene_drugs, real_gene_drugs, max_workers=4, max_papers_per_gene=5,
    )
    assert set(result["gene"].unique()) == {"PTEN"}, "EGFR and TYK2 are well-covered and must not be searched"


def test_build_merged_pool_with_provenance_tracks_pmids(curated_gene_drugs):
    discovery_result = pd.DataFrame([
        {"gene": "PTEN", "candidate_drug": "capivasertib", "pmid": "999", "dgidb_confirmed": True, "abstract_snippet": "..."},
    ])
    merged, prov = provenance.build_merged_pool_with_provenance(discovery_result, curated_gene_drugs)
    assert "capivasertib" in merged["PTEN"]
    assert prov["PTEN"]["capivasertib"] == ["999"]
    assert prov["PTEN"]["alpelisib"] == [], "curated (non-discovered) drug must have empty provenance"


def test_build_merged_pool_with_provenance_rejects_malformed_schema():
    with pytest.raises(ValueError):
        provenance.build_merged_pool_with_provenance(pd.DataFrame(), {"EGFR": "afatinib"})
