"""
tests/test_redundancy_analyzer.py
Unit tests for the RTK/nRTK redundancy analysis module.
"""

import os
import sys
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from modules.vcf_parser import VCFParser
from modules.redundancy_analyzer import RedundancyAnalyzer

SAMPLE_VCF = os.path.join(os.path.dirname(__file__), "..", "data", "sample", "patient_1.vcf")


@pytest.fixture(scope="module")
def analyzer():
    parser = VCFParser(SAMPLE_VCF).parse()
    return RedundancyAnalyzer(parser.rtk_nrtk_variants)


class TestRedundancyAnalyzer:

    def test_affected_genes_populated(self, analyzer):
        assert len(analyzer.affected_genes) > 0

    def test_pathway_map_built(self, analyzer):
        assert len(analyzer._pathway_map) > 0

    def test_redundancy_scores_range(self, analyzer):
        scores = analyzer.compute_redundancy_score()
        for pathway, score in scores.items():
            assert 0.0 <= score <= 1.0, f"{pathway} score {score} out of range"

    def test_mapk_erk_high_redundancy(self, analyzer):
        scores = analyzer.compute_redundancy_score()
        assert "MAPK/ERK" in scores
        assert scores["MAPK/ERK"] > 0.7, "MAPK/ERK should be highly redundant"

    def test_bypass_risk_returns_list(self, analyzer):
        risks = analyzer.bypass_risk_per_target()
        assert isinstance(risks, list)
        assert len(risks) > 0

    def test_bypass_risk_fields(self, analyzer):
        for item in analyzer.bypass_risk_per_target():
            assert "gene" in item
            assert "bypass_risk_score" in item
            assert "recommendation" in item
            assert "bypass_genes" in item
            assert 0.0 <= item["bypass_risk_score"] <= 1.0

    def test_bypass_risk_sorted_descending(self, analyzer):
        risks = analyzer.bypass_risk_per_target()
        scores = [r["bypass_risk_score"] for r in risks]
        assert scores == sorted(scores, reverse=True)

    def test_combination_suggestions_returned(self, analyzer):
        combos = analyzer.combination_therapy_suggestions()
        assert isinstance(combos, list)
        assert len(combos) > 0

    def test_combination_fields(self, analyzer):
        for combo in analyzer.combination_therapy_suggestions():
            assert "target_gene_1" in combo
            assert "target_gene_2" in combo
            assert "shared_pathway" in combo
            assert "rationale" in combo
            assert "evidence_strength" in combo
            assert combo["target_gene_1"] != combo["target_gene_2"]

    def test_full_report_structure(self, analyzer):
        report = analyzer.get_full_report()
        required = {
            "affected_genes",
            "active_pathways",
            "pathway_redundancy_scores",
            "bypass_risk_per_target",
            "combination_therapy_suggestions",
        }
        assert required.issubset(report.keys())

    def test_no_self_bypass(self, analyzer):
        for item in analyzer.bypass_risk_per_target():
            gene = item["gene"]
            bypass_genes = [b["gene"] for b in item["bypass_genes"]]
            assert gene not in bypass_genes, f"{gene} should not list itself as bypass"
