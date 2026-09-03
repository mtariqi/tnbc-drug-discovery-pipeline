"""
diagnose_negative_scores.py

Decomposes DrugPairScore into its real component contributions for all
three FunMap-comparison models, to isolate whether the negative mean
scores observed in the real run come from the redundancy penalty, the
toxicity term, or both -- rather than leaving this unexplained.

Mathematical fact checked here, not just asserted: adjusted_base is
clamped at a minimum of 0.0 (max(base_score - penalty, 0.0)), so the
redundancy penalty alone CANNOT make a pair's contribution negative --
only the toxicity term (-lambda_toxicity * tox) can. The redundancy
penalty can only make a pair more VULNERABLE to a toxicity penalty by
lowering its adjusted_base closer to zero first.

Usage (run from src/scoring/, same directory as the other scripts):
    python3 diagnose_negative_scores.py <path_to_real_funmap.tsv> [base_dir]
"""

import itertools
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))

from funmap_dual_network_scoring import DUAL_NETWORK_MODELS, load_funmap_network, pair_cts_dual_network
from kinase_scoring_pipeline import DRUG_PAIR_WEIGHTS

try:
    from run_pairs_and_triplets import load_real_inputs, _normalize_drug_name
except ImportError:
    raise SystemExit("run_pairs_and_triplets.py must be in the same directory.")

from compare_funmap_vs_string_models import build_toxicity_lookup_from_local_faers


def decompose(
    drugs, target_map, cts, community_map, string_edges, funmap_edges, weights,
    toxicity_lookup, redundancy_penalty_weight=0.15,
):
    """Same two-pass logic as rank_all_pairs_dual_redundancy_penalized(), but
    keeps every intermediate component instead of only the final score."""
    from collections import Counter

    pass_one = []
    winning_pair_counts = Counter()
    for d1, d2 in itertools.combinations(sorted(set(drugs)), 2):
        targets_i = target_map.get(d1, [])
        targets_j = target_map.get(d2, [])
        if not targets_i or not targets_j:
            continue
        best_score, best_pair = float("-inf"), None
        for ki, kj in itertools.product(targets_i, targets_j):
            s = pair_cts_dual_network(ki, kj, cts, community_map, string_edges, funmap_edges, weights)
            if s > best_score:
                best_score, best_pair = s, (ki, kj)
        normalized_pair = tuple(sorted(best_pair)) if best_pair else None
        pass_one.append({"drug_1": d1, "drug_2": d2, "raw_base": best_score, "winning_kinase_pair": normalized_pair})
        if normalized_pair is not None:
            winning_pair_counts[normalized_pair] += 1

    rows = []
    for r in pass_one:
        freq = winning_pair_counts[r["winning_kinase_pair"]] if r["winning_kinase_pair"] else 1
        penalty = redundancy_penalty_weight * max(freq - 1, 0)
        adjusted_base = max(r["raw_base"] - penalty, 0.0)

        d1, d2 = r["drug_1"], r["drug_2"]
        tox = toxicity_lookup.get((d1, d2), toxicity_lookup.get((d2, d1), 0.0))
        tox_contribution = -DRUG_PAIR_WEIGHTS["lambda_toxicity"] * tox
        final_score = adjusted_base + tox_contribution

        rows.append({
            "drug_1": d1, "drug_2": d2,
            "raw_base": r["raw_base"],
            "redundancy_penalty": penalty,
            "adjusted_base": adjusted_base,
            "toxicity_value": tox,
            "toxicity_contribution": tox_contribution,
            "final_score": final_score,
        })
    return pd.DataFrame(rows)


def main():
    if len(sys.argv) < 2:
        raise SystemExit("Usage: python3 diagnose_negative_scores.py <path_to_real_funmap.tsv> [base_dir]")
    funmap_path = sys.argv[1]
    base_dir = sys.argv[2] if len(sys.argv) > 2 else str(Path.home() / "rtk_nrtk_tnbc")

    print("Loading real inputs...")
    cts, community_map, string_edges, target_map = load_real_inputs(base_dir)
    kinase_panel = list(cts.index)
    funmap_edges = load_funmap_network(funmap_path, kinase_panel=kinase_panel)

    drug_list_path = str(Path(base_dir) / "data/raw/drugs/drug_list.txt")
    drug_list = [_normalize_drug_name(d.strip()) for d in open(drug_list_path) if d.strip()]
    scoreable_drugs = [d for d in drug_list if d in target_map and any(k in cts.index for k in target_map[d])]
    max_drugs = 60
    if len(scoreable_drugs) > max_drugs:
        drug_max_cts = {d: max(cts.get(k, 0) for k in target_map[d]) for d in scoreable_drugs}
        scoreable_drugs = sorted(scoreable_drugs, key=lambda d: drug_max_cts[d], reverse=True)[:max_drugs]

    toxicity_lookup = build_toxicity_lookup_from_local_faers(str(Path(base_dir) / "data/processed/faers"), scoreable_drugs)

    print(f"\n{'Model':<25} {'raw_base':>10} {'redund.pen':>11} {'adj.base':>10} {'tox.contrib':>12} {'final':>10}   pct_negative")
    for model_name, weights in DUAL_NETWORK_MODELS.items():
        df = decompose(scoreable_drugs, target_map, cts, community_map, string_edges, funmap_edges, weights, toxicity_lookup)
        pct_negative = (df["final_score"] < 0).mean() * 100
        print(f"{model_name:<25} {df['raw_base'].mean():>10.4f} {df['redundancy_penalty'].mean():>11.4f} "
              f"{df['adjusted_base'].mean():>10.4f} {df['toxicity_contribution'].mean():>12.4f} "
              f"{df['final_score'].mean():>10.4f}   {pct_negative:.1f}%")

        # Direct test of the analytic claim: among pairs with a negative final score,
        # is adjusted_base always > 0 (i.e., toxicity alone flipped it negative),
        # or does adjusted_base itself ever equal 0 (redundancy penalty fully zeroed it)?
        negative_rows = df[df["final_score"] < 0]
        if len(negative_rows) > 0:
            zeroed_by_redundancy = (negative_rows["adjusted_base"] == 0).sum()
            print(f"  -> of {len(negative_rows)} negative-scoring pairs: {zeroed_by_redundancy} had adjusted_base "
                  f"fully zeroed by the redundancy penalty; {len(negative_rows) - zeroed_by_redundancy} had a "
                  f"positive adjusted_base that toxicity alone pushed negative.")
        print(f"  -> mean toxicity value among negative-scoring pairs: {negative_rows['toxicity_value'].mean():.3f}"
              if len(negative_rows) else "  -> no negative-scoring pairs")

    print("\nConclusion check: redundancy_penalty alone cannot cause a negative score (adjusted_base is clamped >= 0).")
    print("If most negative-scoring pairs have a positive adjusted_base, toxicity is confirmed as the necessary")
    print("and sufficient cause; the redundancy penalty's role is only to make pairs MORE vulnerable to it by")
    print("lowering adjusted_base first, not to cause negativity on its own.")


if __name__ == "__main__":
    main()
