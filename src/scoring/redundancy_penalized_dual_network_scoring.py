"""
redundancy_penalized_dual_network_scoring.py

Extends the SAME redundancy-penalty mechanism already adopted for the
STRING-only PairCTS ranking (redundancy_penalized_pair_scoring.py, weight
0.15, documented in the manuscript's Section 3.8) to the dual-network
(STRING/FunMap/Combined) scoring functions in funmap_dual_network_scoring.py.

WHY THIS IS NEEDED: the real three-model comparison run (compare_funmap_vs_
string_models.py) showed the identical ceiling-effect ties under all three
models (one repeated score dominating the top-15 in each) -- the same
underlying max-aggregation problem already confirmed and partially
remediated for the STRING-only case, now confirmed to recur identically
for FunMap-only and Combined scoring, since none of the three currently use
the redundancy-penalized ranking.

MECHANISM (identical to the already-adopted STRING-only version, just
applied to pair_cts_dual_network() instead of pair_cts()): track which
kinase pair wins best_target_pair_cts_dual()'s search for every drug pair;
penalize each pair's score by weight * (frequency - 1), where frequency is
how often that same winning kinase pair recurs across the whole ranking.
A uniquely-winning kinase pair gets frequency=1, penalty=0 -- unchanged.

Default weight matches the adopted manuscript default (0.15), not the
unadopted first-attempt default (0.05) from the original module.
"""

from __future__ import annotations

import itertools
from collections import Counter
from typing import Dict, List, Optional, Tuple

import pandas as pd

from funmap_dual_network_scoring import pair_cts_dual_network, DUAL_NETWORK_MODELS
from kinase_scoring_pipeline import DRUG_PAIR_WEIGHTS


def best_target_pair_cts_dual_with_winner(
    drug_i: str,
    drug_j: str,
    target_map: Dict[str, List[str]],
    cts: pd.Series,
    community_map: Dict[str, int],
    string_edges: Dict[Tuple[str, str], float],
    funmap_edges: Dict[Tuple[str, str], float],
    weights: Dict[str, float],
) -> Tuple[float, Optional[Tuple[str, str]]]:
    """Identical computation to best_target_pair_cts_dual(), but also returns
    WHICH kinase pair won, needed for the redundancy penalty."""
    targets_i = target_map.get(drug_i, [])
    targets_j = target_map.get(drug_j, [])
    if not targets_i or not targets_j:
        return 0.0, None

    best_score = float("-inf")
    best_pair = None
    for ki, kj in itertools.product(targets_i, targets_j):
        s = pair_cts_dual_network(ki, kj, cts, community_map, string_edges, funmap_edges, weights)
        if s > best_score:
            best_score = s
            best_pair = (ki, kj)
    return best_score, best_pair


def rank_all_pairs_dual_redundancy_penalized(
    drugs,
    target_map: Dict[str, List[str]],
    cts: pd.Series,
    community_map: Dict[str, int],
    string_edges: Dict[Tuple[str, str], float],
    funmap_edges: Dict[Tuple[str, str], float],
    weights: Dict[str, float],
    toxicity_lookup: Dict[Tuple[str, str], float],
    synergy_lookup: Dict[Tuple[str, str], float],
    redundancy_penalty_weight: float = 0.15,
) -> pd.DataFrame:
    """Two-pass ranking, identical logic to rank_all_pairs_redundancy_penalized()
    but for any of the three dual-network weight schemes (Model 1/2/3)."""
    pass_one = []
    winning_pair_counts: Counter = Counter()
    for d1, d2 in itertools.combinations(sorted(set(drugs)), 2):
        base_score, winning_pair = best_target_pair_cts_dual_with_winner(
            d1, d2, target_map, cts, community_map, string_edges, funmap_edges, weights
        )
        normalized_pair = tuple(sorted(winning_pair)) if winning_pair else None
        pass_one.append({
            "drug_1": d1, "drug_2": d2, "base_score": base_score, "winning_kinase_pair": normalized_pair
        })
        if normalized_pair is not None:
            winning_pair_counts[normalized_pair] += 1

    rows = []
    for r in pass_one:
        freq = winning_pair_counts[r["winning_kinase_pair"]] if r["winning_kinase_pair"] else 1
        penalty = redundancy_penalty_weight * max(freq - 1, 0)
        adjusted_base = max(r["base_score"] - penalty, 0.0)

        d1, d2 = r["drug_1"], r["drug_2"]
        tox = toxicity_lookup.get((d1, d2), toxicity_lookup.get((d2, d1), 0.0))
        syn = synergy_lookup.get((d1, d2), synergy_lookup.get((d2, d1), 0.0))
        score = adjusted_base - DRUG_PAIR_WEIGHTS["lambda_toxicity"] * tox + DRUG_PAIR_WEIGHTS["mu_synergy"] * syn

        rows.append({
            "drug_1": d1, "drug_2": d2, "DrugPairScore": score,
            "winning_kinase_pair": r["winning_kinase_pair"], "winning_pair_frequency": freq,
        })

    return pd.DataFrame(rows).sort_values("DrugPairScore", ascending=False).reset_index(drop=True)


# =====================================================================
# SMOKE TEST -- same synthetic-case discipline as the original module:
# confirms the mechanism actually de-ties a case designed to need it,
# and confirms a uniquely-winning pair is genuinely untouched.
# =====================================================================

def _run_smoke_test():
    cts = pd.Series({"ERBB2": 0.9, "EPHA5": 0.6, "EGFR": 0.5, "SRC": 0.3, "ABL1": 0.2})
    community_map = {"ERBB2": 0, "EPHA5": 1, "EGFR": 0, "SRC": 1, "ABL1": 1}
    string_edges = {}
    funmap_edges = {}
    weights = DUAL_NETWORK_MODELS["model2_funmap_only"]

    # Three drugs whose targets all resolve to the same winning kinase pair
    # (ERBB2, EPHA5) against a shared partner -- the same synthetic pattern
    # already used to validate the STRING-only version.
    target_map = {
        "drugA": ["ERBB2"], "drugB": ["ERBB2"], "drugC": ["ERBB2"],
        "partner": ["EPHA5"],
        "uniqueDrug": ["SRC"], "uniquePartner": ["ABL1"],
    }
    drugs = ["drugA", "drugB", "drugC", "partner", "uniqueDrug", "uniquePartner"]
    toxicity_lookup, synergy_lookup = {}, {}

    print("=== Testing redundancy penalty de-ties the repeated (ERBB2, EPHA5) winner ===")
    ranked = rank_all_pairs_dual_redundancy_penalized(
        drugs, target_map, cts, community_map, string_edges, funmap_edges, weights,
        toxicity_lookup, synergy_lookup, redundancy_penalty_weight=0.15,
    )
    erbb2_epha5_rows = ranked[ranked["winning_kinase_pair"] == ("EPHA5", "ERBB2")]
    assert len(erbb2_epha5_rows) == 3, f"expected 3 pairs winning on (EPHA5,ERBB2), got {len(erbb2_epha5_rows)}"
    assert erbb2_epha5_rows["winning_pair_frequency"].iloc[0] == 3
    scores = erbb2_epha5_rows["DrugPairScore"].tolist()
    # All three should be penalized identically (same frequency), so still tied
    # with each other, but strictly lower than their unpenalized base score.
    unpenalized_base = 0.35 * cts["ERBB2"] + 0.35 * cts["EPHA5"] + 0.20 * 1.0  # cross-community, no edges
    for s in scores:
        assert s < unpenalized_base - 1e-9, f"expected penalized score below {unpenalized_base}, got {s}"
    print(f"PASSED: all 3 repeated-winner pairs penalized identically "
          f"(base={unpenalized_base:.4f}, penalized={scores[0]:.4f}).\n")

    print("=== Testing a uniquely-winning kinase pair is NOT penalized ===")
    unique_row = ranked[ranked["winning_kinase_pair"] == ("ABL1", "SRC")]
    assert len(unique_row) == 1
    assert unique_row["winning_pair_frequency"].iloc[0] == 1
    expected_unique_score = 0.35 * cts["SRC"] + 0.35 * cts["ABL1"] + 0.20 * 0.0  # same community (both = 1)
    actual = unique_row["DrugPairScore"].iloc[0]
    assert abs(actual - expected_unique_score) < 1e-9, f"expected {expected_unique_score}, got {actual}"
    print(f"PASSED: uniquely-winning pair scored {actual:.4f}, exactly matching its unpenalized base score "
          f"(frequency=1 -> penalty=0).\n")

    print("=" * 70)
    print("ALL SMOKE TESTS PASSED. Ready to apply to the real three-model comparison.")


if __name__ == "__main__":
    _run_smoke_test()
