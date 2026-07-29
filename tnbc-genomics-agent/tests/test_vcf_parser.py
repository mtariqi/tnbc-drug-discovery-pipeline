"""
tests/test_vcf_parser.py
Unit tests for the VCF parsing module.
"""

import os
import sys
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from modules.vcf_parser import VCFParser, Variant

SAMPLE_VCF = os.path.join(os.path.dirname(__file__), "..", "data", "sample", "patient_1.vcf")


class TestVCFParser:

    def test_parse_returns_self(self):
        parser = VCFParser(SAMPLE_VCF)
        result = parser.parse()
        assert result is parser

    def test_total_variants_parsed(self):
        parser = VCFParser(SAMPLE_VCF).parse()
        assert parser.get_summary()["total_variants"] == 13

    def test_rtk_nrtk_variants_detected(self):
        parser = VCFParser(SAMPLE_VCF).parse()
        assert parser.get_summary()["rtk_nrtk_variants"] == 10

    def test_known_genes_detected(self):
        parser = VCFParser(SAMPLE_VCF).parse()
        genes = set(parser.get_summary()["genes_detected"])
        expected = {"EGFR", "MET", "SRC", "ABL1", "FGFR1", "FGFR2", "FGFR3", "ERBB2", "ALK"}
        assert expected == genes

    def test_variant_fields_populated(self):
        parser = VCFParser(SAMPLE_VCF).parse()
        for v in parser.variants:
            assert isinstance(v, Variant)
            assert v.chrom
            assert v.pos > 0
            assert v.gene
            assert 0.0 <= v.allele_freq <= 1.0
            assert v.depth >= 0

    def test_rtk_nrtk_flag_set_correctly(self):
        parser = VCFParser(SAMPLE_VCF).parse()
        for v in parser.rtk_nrtk_variants:
            assert v.is_rtk_nrtk is True
            assert v.kinase_info is not None

    def test_severity_scores_valid(self):
        parser = VCFParser(SAMPLE_VCF).parse()
        for v in parser.variants:
            assert 0.0 <= v.severity_score <= 1.0

    def test_high_severity_variants_in_summary(self):
        parser = VCFParser(SAMPLE_VCF).parse()
        high = parser.get_summary()["high_severity_variants"]
        assert len(high) > 0
        for v in high:
            assert v["severity_score"] >= 0.7

    def test_to_dict_keys(self):
        parser = VCFParser(SAMPLE_VCF).parse()
        v = parser.variants[0]
        d = v.to_dict()
        required_keys = {"chrom", "pos", "gene", "effect", "allele_freq",
                         "depth", "severity_score", "is_rtk_nrtk"}
        assert required_keys.issubset(d.keys())

    def test_missing_file_raises(self):
        with pytest.raises(FileNotFoundError):
            VCFParser("nonexistent.vcf").parse()
