"""
Compare MDCOE Results: Curated vs. Real DGIdb-Sourced Drug Candidates
========================================================================

Runs the SAME patient/gene profile through MDCOE twice -- once using
mdcoe.py's hand-curated GENE_DRUGS (the path that already produced your
validated afatinib+alpelisib+trastuzumab, HCOS=0.450 result), and once
using REAL_GENE_DRUGS built from live DGIdb data (build_real_gene_drugs.py)
-- and reports both side by side. Does NOT replace or modify the curated
path; your already-validated result remains exactly reproducible.

ONE IMPORTANT DIFFERENCE TO EXPECT, STATED UP FRONT:
    TP53 is handled specially in the curated path (resolve_tp53_drugs()
    maps a TP53 alteration to WEE1/BCL-2 inhibitors based on the biology
    of NMD vs. missense mutations -- there is no direct, validated TP53
    inhibitor for most alteration types). The real-DGIdb path has no such
    biological special-casing -- it will return whatever DGIdb happens to
    have annotated for the literal gene symbol "TP53" (which may include
    noisy chemotherapy-agent annotations, or real reactivator compounds
    like eprenetapopt, or nothing at all). This is a genuine, expected
    difference between the two approaches, not a bug -- don't force them
    to match.
"""

from __future__ import annotations

from typing import Dict, List

from mdcoe import GENE_DRUGS, resolve_tp53_drugs, DrugGraph, SynergyNet, HCOS, MDCOE
from build_real_gene_drugs import build_real_gene_drugs


def build_curated_candidates(genes: List[str], tp53_alteration_type: str = "Nonsense_Mutation") -> List[str]:
    """The existing, already-validated candidate-drug derivation path."""
    drugs = []
    for g in genes:
        if g == "TP53":
            drugs.extend(resolve_tp53_drugs(tp53_alteration_type))
        else:
            drugs.extend(GENE_DRUGS.get(g, []))
    return sorted(set(drugs))


def build_real_candidates(
    genes: List[str],
    dgidb_interactions_path: str,
    all_kinase_panel: List[str],
    max_drugs_per_gene: int = 3,
) -> List[str]:
    """The new, real-DGIdb-sourced candidate-drug derivation path."""
    real_gene_drugs = build_real_gene_drugs(
        dgidb_interactions_path, genes, max_drugs_per_gene=max_drugs_per_gene,
        all_kinase_panel=all_kinase_panel,
    )
    drugs = []
    for g in genes:
        drugs.extend(real_gene_drugs.get(g, []))
    return sorted(set(drugs))


def compare(
    genes: List[str],
    dgidb_interactions_path: str,
    all_kinase_panel: List[str],
    tp53_alteration_type: str = "Nonsense_Mutation",
    top_k: int = 10,
):
    curated_drugs = build_curated_candidates(genes, tp53_alteration_type)
    real_drugs = build_real_candidates(genes, dgidb_interactions_path, all_kinase_panel)

    print(f"Genes: {genes}")
    print(f"Curated candidates ({len(curated_drugs)}): {curated_drugs}")
    print(f"Real DGIdb-sourced candidates ({len(real_drugs)}): {real_drugs}")

    only_curated = set(curated_drugs) - set(real_drugs)
    only_real = set(real_drugs) - set(curated_drugs)
    if only_curated:
        print(f"\nIn curated but NOT in real-DGIdb candidates (DGIdb may not have this interaction annotated, "
              f"or it's the TP53-special-case biology -- see module docstring): {sorted(only_curated)}")
    if only_real:
        print(f"\nIn real-DGIdb but NOT in curated candidates (potentially new candidates worth reviewing): {sorted(only_real)}")

    net = SynergyNet()

    print(f"\n{'='*70}\nCURATED path -- top {top_k}\n{'='*70}")
    curated_results = MDCOE(DrugGraph(curated_drugs), net, HCOS, beam_width=50, max_depth=5, top_k=top_k)
    for i, (regimen, score) in enumerate(curated_results, 1):
        print(f"  {i:2d}. HCOS={score:.3f}  ({len(regimen)}-drug)  {' + '.join(sorted(regimen))}")

    print(f"\n{'='*70}\nREAL DGIdb-sourced path -- top {top_k}\n{'='*70}")
    if len(real_drugs) < 2:
        print("  Fewer than 2 real candidate drugs -- not enough to search. "
              "This can happen if DGIdb's annotations for these specific genes "
              "were mostly filtered out as noise, or max_drugs_per_gene is too low.")
        real_results = []
    else:
        real_results = MDCOE(DrugGraph(real_drugs), net, HCOS, beam_width=50, max_depth=5, top_k=top_k)
        for i, (regimen, score) in enumerate(real_results, 1):
            print(f"  {i:2d}. HCOS={score:.3f}  ({len(regimen)}-drug)  {' + '.join(sorted(regimen))}")

    return curated_results, real_results


if __name__ == "__main__":
    import sys
    from pathlib import Path

    base_dir = Path(sys.argv[1] if len(sys.argv) > 1 else str(Path.home() / "rtk_nrtk_tnbc"))
    genes = ["EGFR", "PTEN", "FLT1", "TP53", "ERBB2", "TYK2"]  # TCGA-AO-A128
    kinase_panel = [k.strip() for k in open(base_dir / "data/raw/kinases/kinase_90_list.txt")]
    dgidb_path = str(base_dir / "data/processed/dgidb/dgidb_interactions.tsv")

    compare(genes, dgidb_path, kinase_panel)
