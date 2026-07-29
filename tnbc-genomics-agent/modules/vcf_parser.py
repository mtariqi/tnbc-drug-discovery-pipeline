"""
vcf_parser.py
Parses VCF files and extracts RTK/nRTK-relevant variants for TNBC analysis.
"""

import re
from dataclasses import dataclass, field
from typing import List, Dict, Optional
from modules.kinase_db import KINASE_DATABASE, EFFECT_SEVERITY


@dataclass
class Variant:
    chrom: str
    pos: int
    variant_id: str
    ref: str
    alt: str
    qual: str
    filter_status: str
    depth: int
    allele_freq: float
    gene: str
    effect: str
    severity_score: float = 0.0
    is_rtk_nrtk: bool = False
    kinase_info: Optional[Dict] = field(default=None)

    def to_dict(self) -> Dict:
        return {
            "chrom": self.chrom,
            "pos": self.pos,
            "variant_id": self.variant_id,
            "ref": self.ref,
            "alt": self.alt,
            "qual": self.qual,
            "filter_status": self.filter_status,
            "depth": self.depth,
            "allele_freq": self.allele_freq,
            "gene": self.gene,
            "effect": self.effect,
            "severity_score": self.severity_score,
            "is_rtk_nrtk": self.is_rtk_nrtk,
            "kinase_type": self.kinase_info.get("type", "N/A") if self.kinase_info else "N/A",
            "kinase_family": self.kinase_info.get("family", "N/A") if self.kinase_info else "N/A",
            "tnbc_relevance": self.kinase_info.get("tnbc_relevance", "N/A") if self.kinase_info else "N/A",
        }


class VCFParser:
    """
    Parses a standard VCF file and identifies RTK/nRTK variants relevant to TNBC.
    Supports GENE and EFFECT fields in the INFO column.
    """

    def __init__(self, vcf_path: str):
        self.vcf_path = vcf_path
        self.variants: List[Variant] = []
        self.rtk_nrtk_variants: List[Variant] = []
        self.metadata: Dict = {}
        self._known_genes = set(KINASE_DATABASE.keys())

    def parse(self) -> "VCFParser":
        """Parse the VCF file and populate variant lists."""
        with open(self.vcf_path, "r") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith("##"):
                    self._parse_meta(line)
                    continue
                if line.startswith("#CHROM"):
                    self.metadata["samples"] = line.split("\t")[9:]
                    continue
                variant = self._parse_record(line)
                if variant:
                    self.variants.append(variant)
                    if variant.is_rtk_nrtk:
                        self.rtk_nrtk_variants.append(variant)
        return self

    # ── Private helpers ──────────────────────────────────────────────────────

    def _parse_meta(self, line: str):
        """Store key VCF metadata."""
        if "fileformat" in line:
            self.metadata["fileformat"] = line.split("=", 1)[1]
        elif "reference" in line.lower():
            self.metadata["reference"] = line.split("=", 1)[1]

    def _parse_record(self, line: str) -> Optional[Variant]:
        """Parse a single VCF data line into a Variant object."""
        cols = line.split("\t")
        if len(cols) < 8:
            return None

        chrom, pos, vid, ref, alt, qual, filt, info_str = cols[:8]
        info = self._parse_info(info_str)

        gene   = info.get("GENE", self._guess_gene_from_pos(chrom, int(pos)))
        effect = info.get("EFFECT", "unknown")
        depth  = int(info.get("DP", 0))
        af     = float(info.get("AF", 0.0))

        severity = EFFECT_SEVERITY.get(effect, 0.0)
        is_kinase = gene in self._known_genes
        kinase_info = KINASE_DATABASE.get(gene) if is_kinase else None

        return Variant(
            chrom=chrom,
            pos=int(pos),
            variant_id=vid,
            ref=ref,
            alt=alt,
            qual=qual,
            filter_status=filt,
            depth=depth,
            allele_freq=af,
            gene=gene,
            effect=effect,
            severity_score=severity,
            is_rtk_nrtk=is_kinase,
            kinase_info=kinase_info,
        )

    @staticmethod
    def _parse_info(info_str: str) -> Dict:
        """Convert semicolon-delimited INFO field into a dict."""
        result: Dict = {}
        for token in info_str.split(";"):
            if "=" in token:
                k, v = token.split("=", 1)
                result[k] = v
            else:
                result[token] = True
        return result

    @staticmethod
    def _guess_gene_from_pos(chrom: str, pos: int) -> str:
        """Fallback: approximate gene name from chromosomal coordinates."""
        # Minimal locus map for demo purposes
        locus_map = [
            ("7",  55_000_000,  55_300_000, "EGFR"),
            ("17", 37_800_000,  37_900_000, "ERBB2"),
            ("7",  116_300_000, 116_500_000,"MET"),
            ("9",  107_500_000, 107_600_000,"ABL1"),
            ("20", 37_300_000,  37_500_000, "SRC"),
        ]
        for c, start, end, gene in locus_map:
            if chrom == c and start <= pos <= end:
                return gene
        return "UNKNOWN"

    # ── Accessors ────────────────────────────────────────────────────────────

    def get_summary(self) -> Dict:
        """Return a high-level parse summary."""
        return {
            "total_variants": len(self.variants),
            "rtk_nrtk_variants": len(self.rtk_nrtk_variants),
            "genes_detected": list({v.gene for v in self.rtk_nrtk_variants}),
            "high_severity_variants": [
                v.to_dict() for v in self.rtk_nrtk_variants if v.severity_score >= 0.7
            ],
        }
