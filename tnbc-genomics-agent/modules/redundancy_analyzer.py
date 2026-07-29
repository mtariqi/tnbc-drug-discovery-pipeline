"""
redundancy_analyzer.py
Detects RTK/nRTK pathway redundancy in a patient's variant profile.
Scores the risk that blocking one kinase will be bypassed by another.
"""

from typing import List, Dict, Tuple
from collections import defaultdict
from modules.kinase_db import KINASE_DATABASE, PATHWAY_REDUNDANCY_MAP
from modules.vcf_parser import Variant


class RedundancyAnalyzer:
    """
    Given a list of RTK/nRTK variants from a patient, computes:
      - Active pathways
      - Pathway overlap (redundancy) score
      - Bypass risk per gene target
      - Combination therapy suggestions
    """

    def __init__(self, rtk_nrtk_variants: List[Variant]):
        self.variants = rtk_nrtk_variants
        self.affected_genes = {v.gene for v in rtk_nrtk_variants}
        self._pathway_map: Dict[str, List[str]] = defaultdict(list)
        self._build_pathway_map()

    # ── Core analysis ────────────────────────────────────────────────────────

    def _build_pathway_map(self):
        """Map each affected pathway to the genes feeding into it."""
        for gene in self.affected_genes:
            info = KINASE_DATABASE.get(gene)
            if not info:
                continue
            for pathway in info.get("pathways", []):
                self._pathway_map[pathway].append(gene)

    def compute_redundancy_score(self) -> Dict[str, float]:
        """
        Returns a per-pathway redundancy score (0–1).
        Score reflects how many independent genes feed that pathway.
        Higher = more redundant = harder to block.
        """
        scores: Dict[str, float] = {}
        for pathway, genes in self._pathway_map.items():
            n = len(genes)
            # Sigmoid-like scaling: 1 gene → 0.2, 2 → 0.5, 3 → 0.75, 4+ → >0.85
            score = 1 - (1 / (1 + 0.8 * (n - 1)))
            scores[pathway] = round(score, 3)
        return scores

    def bypass_risk_per_target(self) -> List[Dict]:
        """
        For each affected RTK/nRTK gene, estimate the bypass risk:
        how many other affected genes could rescue the same pathway if
        this gene were inhibited.
        """
        results = []
        for gene in sorted(self.affected_genes):
            info = KINASE_DATABASE.get(gene)
            if not info:
                continue

            gene_pathways = set(info.get("pathways", []))
            # Find OTHER affected genes sharing at least one pathway
            bypass_genes = []
            for other_gene in self.affected_genes:
                if other_gene == gene:
                    continue
                other_info = KINASE_DATABASE.get(other_gene)
                if not other_info:
                    continue
                shared = gene_pathways & set(other_info.get("pathways", []))
                if shared:
                    bypass_genes.append({
                        "gene": other_gene,
                        "shared_pathways": list(shared),
                        "type": other_info["type"],
                    })

            # Also check known resistance mechanisms
            resistance_hits = [
                r for r in info.get("resistance_mechanisms", [])
                if any(r.lower().find(ag.lower()) != -1 for ag in self.affected_genes)
            ]

            n_bypass = len(bypass_genes)
            risk_score = min(1.0, round(0.15 * n_bypass + 0.25 * len(resistance_hits), 3))

            results.append({
                "gene": gene,
                "type": info["type"],
                "family": info["family"],
                "tnbc_relevance": info["tnbc_relevance"],
                "pathways": list(gene_pathways),
                "bypass_genes": bypass_genes,
                "resistance_mechanisms_present": resistance_hits,
                "bypass_risk_score": risk_score,
                "recommendation": self._risk_label(risk_score),
            })

        return sorted(results, key=lambda x: x["bypass_risk_score"], reverse=True)

    def combination_therapy_suggestions(self) -> List[Dict]:
        """
        Suggests rational drug combinations that could co-inhibit redundant pathways.
        Groups genes by shared pathway and proposes paired inhibitor strategies.
        """
        suggestions = []
        seen_pairs = set()

        for pathway, genes in self._pathway_map.items():
            if len(genes) < 2:
                continue

            for i in range(len(genes)):
                for j in range(i + 1, len(genes)):
                    gA, gB = genes[i], genes[j]
                    pair_key = tuple(sorted([gA, gB]))
                    if pair_key in seen_pairs:
                        continue
                    seen_pairs.add(pair_key)

                    drugs_A = [d["drug"] for d in KINASE_DATABASE[gA].get("known_inhibitors", [])]
                    drugs_B = [d["drug"] for d in KINASE_DATABASE[gB].get("known_inhibitors", [])]

                    # Prefer approved drugs
                    approved_A = [d["drug"] for d in KINASE_DATABASE[gA].get("known_inhibitors", []) if d.get("fda_approved")]
                    approved_B = [d["drug"] for d in KINASE_DATABASE[gB].get("known_inhibitors", []) if d.get("fda_approved")]

                    suggestions.append({
                        "target_gene_1": gA,
                        "target_gene_2": gB,
                        "shared_pathway": pathway,
                        "rationale": (
                            f"Both {gA} and {gB} activate {pathway}. "
                            f"Inhibiting only {gA} risks compensation via {gB}."
                        ),
                        "suggested_drugs_gene1": approved_A[:2] or drugs_A[:2],
                        "suggested_drugs_gene2": approved_B[:2] or drugs_B[:2],
                        "evidence_strength": self._evidence_strength(gA, gB),
                    })

        return sorted(suggestions, key=lambda x: x["evidence_strength"], reverse=True)

    def get_full_report(self) -> Dict:
        """Aggregate all redundancy analyses into one report dict."""
        return {
            "affected_genes": sorted(self.affected_genes),
            "active_pathways": dict(self._pathway_map),
            "pathway_redundancy_scores": self.compute_redundancy_score(),
            "bypass_risk_per_target": self.bypass_risk_per_target(),
            "combination_therapy_suggestions": self.combination_therapy_suggestions(),
        }

    # ── Helpers ──────────────────────────────────────────────────────────────

    @staticmethod
    def _risk_label(score: float) -> str:
        if score >= 0.6:
            return "HIGH BYPASS RISK — combination therapy strongly recommended"
        elif score >= 0.3:
            return "MODERATE BYPASS RISK — consider combination or sequential therapy"
        else:
            return "LOW BYPASS RISK — monotherapy may be effective"

    @staticmethod
    def _evidence_strength(gA: str, gB: str) -> str:
        """
        Rough evidence tier based on TNBC relevance of both genes.
        """
        rel_map = {"high": 3, "moderate": 2, "low": 1}
        rA = rel_map.get(KINASE_DATABASE.get(gA, {}).get("tnbc_relevance", "low"), 1)
        rB = rel_map.get(KINASE_DATABASE.get(gB, {}).get("tnbc_relevance", "low"), 1)
        score = rA + rB
        if score >= 5:
            return "strong"
        elif score >= 3:
            return "moderate"
        else:
            return "limited"
