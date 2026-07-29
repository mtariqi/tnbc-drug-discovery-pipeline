"""
tests/test_kinase_db.py
Unit tests for the curated kinase database.
"""

import os
import sys
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from modules.kinase_db import KINASE_DATABASE, PATHWAY_REDUNDANCY_MAP, EFFECT_SEVERITY


class TestKinaseDatabase:

    def test_core_tnbc_genes_present(self):
        required = {"EGFR", "MET", "FGFR1", "FGFR2", "FGFR3", "ERBB2", "SRC", "ABL1", "ALK"}
        assert required.issubset(KINASE_DATABASE.keys())

    def test_each_gene_has_required_fields(self):
        required_fields = {
            "full_name", "type", "family", "chromosome",
            "tnbc_relevance", "pathways", "known_inhibitors",
            "resistance_mechanisms", "redundant_pathways",
        }
        for gene, info in KINASE_DATABASE.items():
            missing = required_fields - info.keys()
            assert not missing, f"{gene} missing fields: {missing}"

    def test_type_values_valid(self):
        valid_types = {"RTK", "nRTK"}
        for gene, info in KINASE_DATABASE.items():
            assert info["type"] in valid_types, f"{gene} has invalid type: {info['type']}"

    def test_tnbc_relevance_values_valid(self):
        valid = {"high", "moderate", "low"}
        for gene, info in KINASE_DATABASE.items():
            assert info["tnbc_relevance"] in valid

    def test_inhibitors_have_required_fields(self):
        for gene, info in KINASE_DATABASE.items():
            for inh in info["known_inhibitors"]:
                assert "drug" in inh
                assert "fda_approved" in inh
                assert isinstance(inh["fda_approved"], bool)

    def test_overexpression_freq_in_range(self):
        for gene, info in KINASE_DATABASE.items():
            freq = info.get("overexpression_freq_tnbc", 0)
            assert 0.0 <= freq <= 1.0, f"{gene} freq {freq} out of range"

    def test_pathways_nonempty(self):
        for gene, info in KINASE_DATABASE.items():
            assert len(info["pathways"]) > 0, f"{gene} has no pathways"

    def test_pathway_redundancy_map_nonempty(self):
        assert len(PATHWAY_REDUNDANCY_MAP) > 0

    def test_effect_severity_scores_in_range(self):
        for effect, score in EFFECT_SEVERITY.items():
            assert 0.0 <= score <= 1.0, f"{effect} severity {score} out of range"

    def test_stop_gained_highest_severity(self):
        assert EFFECT_SEVERITY["stop_gained"] >= EFFECT_SEVERITY["missense_variant"]
        assert EFFECT_SEVERITY["stop_gained"] >= EFFECT_SEVERITY["synonymous_variant"]
