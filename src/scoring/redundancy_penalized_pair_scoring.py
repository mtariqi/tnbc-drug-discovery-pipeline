"""
redundancy_penalized_pair_scoring.py

Second attempt at the confirmed PairCTS ceiling-effect problem (see
docs/limitations.md). The first attempt (primary_target_pair_scoring.py) was
tested against real data and made things WORSE (unique scores dropped from
430 to 251) -- documented as a rejected negative result, not deleted.

DIFFERENT APPROACH THIS TIME: rather than changing HOW any single pair's
score is computed (which is what failed last time), this leaves
best_target_pair_cts()'s search completely untouched, and instead tracks
WHICH kinase pair wins that search for every drug pair in the ranking. If
the same winning kinase pair recurs across many different drug pairs (the
actual confirmed problem -- e.g. (ERBB2, EPHA5) winning for both
dacomitinib+hesperadin and crizotinib+hesperadin), each of those pairs'
scores is penalized in proportion to how often that same kinase pair
recurred elsewhere in the ranking. A pair whose winning kinase combination
is unique in the ranking gets NO penalty at all.

HONEST CAVEAT, learned from the first attempt's failure: this is tested
carefully on synthetic data below to confirm the mechanism does what it's
supposed to, but whether it meaningfully helps on the REAL 116-drug dataset
is an empirical question, not a certainty -- exactly the mistake made last
time was asserting the fix would help before checking. Run
compare_redundancy_penalty.py (companion script) against real data before
concluding anything.
"""

from __future__ import annotations

import itertools
from collections import Counter
from typing import Dict, List, Optional, Tuple

import pandas as pd

from kinase_scoring_pipeline import pair_cts, DRUG_PAIR_WEIGHTS


def best_target_pair_cts_with_winner(
    drug_i: str,
    drug_j: str,
    target_map: Dict[str, List[str]],
    cts: pd.Series,
    community_map: Dict[str, int],
    crosstalk_edges: Dict[Tuple[str, str], float],
) -> Tuple[float, Optional[Tuple[str, str]]]:
    """
    Identical computation to kinase_scoring_pipeline.py's best_target_pair_cts()
    -- same score, same tie-breaking (first-encountered max in
    itertools.product(targets_i, targets_j) order, matching Python's max() on
    a generator) -- but also returns WHICH kinase pair won, since that
    information is needed for the redundancy penalty and the original
    function doesn't expose it.
    """
    targets_i = target_map.get(drug_i, [])
    targets_j = target_map.get(drug_j, [])
    if not targets_i or not targets_j:
        return 0.0, None

    best_score = float("-inf")
    best_pair = None
    for ki, kj in itertools.product(targets_i, targets_j):
        s = pair_cts(ki, kj, cts, community_map, crosstalk_edges)
        if s > best_score:
            best_score = s
            best_pair = (ki, kj)
    return best_score, best_pair


def rank_all_pairs_redundancy_penalized(
    drugs,
    target_map: Dict[str, List[str]],
    cts: pd.Series,
    community_map: Dict[str, int],
    crosstalk_edges: Dict[Tuple[str, str], float],
    toxicity_lookup: Dict[Tuple[str, str], float],
    synergy_lookup: Dict[Tuple[str, str], float],
    redundancy_penalty_weight: float = 0.05,
) -> pd.DataFrame:
    """
    Two-pass ranking: first computes every pair's base score AND records which
    kinase pair won (via best_target_pair_cts_with_winner(), identical score
    to the original). Second pass counts how often each winning kinase pair
    recurred across the WHOLE ranking, and subtracts
    redundancy_penalty_weight * (frequency - 1) from each pair's base score --
    a pair whose winning kinase combination is unique gets frequency=1, so
    penalty=0, unchanged from the original. Toxicity/synergy terms are then
    applied exactly as drug_pair_score() does.
    """
    pass_one = []
    winning_pair_counts: Counter = Counter()
    for d1, d2 in itertools.combinations(sorted(set(drugs)), 2):
        base_score, winning_pair = best_target_pair_cts_with_winner(
            d1, d2, target_map, cts, community_map, crosstalk_edges
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
# SMOKE TEST -- verifies the mechanism does what it's supposed to on a
# synthetic case designed to unambiguously demonstrate it, not just run
# without crashing (the lesson from the first attempt's failure)
# =====================================================================

def _run_smoke_test():
    from kinase_scoring_pipeline import best_target_pair_cts, rank_all_pairs

    cts = pd.Series({"ERBB2": 0.9, "EPHA5": 0.6, "EGFR": 0.5, "SRC": 0.3, "ABL1": 0.2})
    community_map = {"ERBB2": 0, "EPHA5": 1, "EGFR": 0, "SRC": 1, "ABL1": 1}
    crosstalk_edges = {}

    # Three drugs that ALL independently resolve to the same winning kinase pair
    # (ERBB2, EPHA5) when paired with a common partner -- mirrors the real
    # confirmed finding (dacomitinib+hesperadin and crizotinib+hesperadin both
    # resolving to (ERBB2, EPHA5)). Plus one drug pair with a genuinely unique
    # winning kinase pair, which should NOT be penalized.
    target_map = {
        "partner_hub": ["EPHA5"],  # plays the role of "hesperadin" in this synthetic case
        "drug_a": ["ERBB2"],
        "drug_b": ["ERBB2"],
        "drug_c": ["ERBB2"],
        "unique_x": ["SRC"],
        "unique_y": ["ABL1"],
    }
    toxicity_lookup, synergy_lookup = {}, {}
    drugs = list(target_map.keys())

    print("=== Testing best_target_pair_cts_with_winner() matches the original's score exactly ===")
    for d1, d2 in itertools.combinations(drugs, 2):
        original_score = best_target_pair_cts(d1, d2, target_map, cts, community_map, crosstalk_edges)
        new_score, winner = best_target_pair_cts_with_winner(d1, d2, target_map, cts, community_map, crosstalk_edges)
        assert abs(original_score - new_score) < 1e-12, \
            f"{d1}+{d2}: original={original_score}, with_winner={new_score} -- scores must match exactly"
    print("PASSED: identical scores to the original for every pair -- the winner-tracking wrapper "
          "doesn't change any base computation, only adds bookkeeping.\n")

    print("=== Testing three drugs sharing the same winning kinase pair are correctly identified ===")
    _, w1 = best_target_pair_cts_with_winner("drug_a", "partner_hub", target_map, cts, community_map, crosstalk_edges)
    _, w2 = best_target_pair_cts_with_winner("drug_b", "partner_hub", target_map, cts, community_map, crosstalk_edges)
    _, w3 = best_target_pair_cts_with_winner("drug_c", "partner_hub", target_map, cts, community_map, crosstalk_edges)
    assert tuple(sorted(w1)) == tuple(sorted(w2)) == tuple(sorted(w3)) == ("EPHA5", "ERBB2")
    print(f"PASSED: drug_a, drug_b, drug_c each independently resolve to the same winning kinase "
          f"pair {tuple(sorted(w1))} when paired with partner_hub, exactly mirroring the real "
          f"confirmed finding.\n")

    print("=== Testing the redundancy penalty actually reduces the repeated pairs' scores ===")
    ranked = rank_all_pairs_redundancy_penalized(
        drugs, target_map, cts, community_map, crosstalk_edges, toxicity_lookup, synergy_lookup,
        redundancy_penalty_weight=0.05,
    )
    print(ranked.to_string(index=False))

    repeated_rows = ranked[ranked["drug_2"] == "partner_hub"]
    repeated_rows = repeated_rows[repeated_rows["drug_1"].isin(["drug_a", "drug_b", "drug_c"])]
    assert (repeated_rows["winning_pair_frequency"] == 3).all(), "all three should show frequency=3"
    original_unpenalized_score = best_target_pair_cts("drug_a", "partner_hub", target_map, cts, community_map, crosstalk_edges)
    assert (repeated_rows["DrugPairScore"] < original_unpenalized_score - 1e-9).all(), \
        "all three repeated pairs should score LOWER than their unpenalized original score"
    expected_penalty = 0.05 * (3 - 1)  # frequency=3, weight=0.05
    assert abs(repeated_rows["DrugPairScore"].iloc[0] - (original_unpenalized_score - expected_penalty)) < 1e-9
    print(f"PASSED: all three drugs sharing the (EPHA5, ERBB2) winning pair are correctly penalized "
          f"down from their original {original_unpenalized_score:.4f} by exactly {expected_penalty:.4f} "
          f"(frequency=3, weight=0.05).\n")

    print("=== Testing a pair with a UNIQUE winning kinase pair is NOT penalized ===")
    unique_row = ranked[(ranked["drug_1"] == "unique_x") & (ranked["drug_2"] == "unique_y")]
    assert len(unique_row) == 1
    assert unique_row["winning_pair_frequency"].iloc[0] == 1
    original_unique_score = best_target_pair_cts("unique_x", "unique_y", target_map, cts, community_map, crosstalk_edges)
    assert abs(unique_row["DrugPairScore"].iloc[0] - original_unique_score) < 1e-9, \
        "a pair with a unique (frequency=1) winning kinase pair should be completely unpenalized"
    print(f"PASSED: unique_x+unique_y (winning pair frequency=1) scores exactly {original_unique_score:.4f}, "
          f"identical to its unpenalized original -- confirming genuinely unique pairs are left alone.\n")

    print("=== Testing rank-order can flip in favor of a unique match over a redundant one ===")
    # A fair test of "does this improve ranking usefulness" isn't "more unique score values"
    # (that metric can legitimately stay flat if some drugs are truly interchangeable by
    # construction, as drug_a/b/c are here -- all three target ONLY ERBB2, so no scoring
    # method could ever tell them apart, redundancy-penalized or not). The real, defensible
    # claim is narrower: a pair whose winning kinase-pair is common should be able to drop
    # BELOW a pair whose winning kinase-pair is unique, even if the unique pair's raw score
    # was originally lower -- letting genuine novelty surface past a repeated ceiling.
    redundant_pair_score = ranked[
        (ranked["drug_1"] == "drug_a") & (ranked["drug_2"] == "partner_hub")
    ]["DrugPairScore"].iloc[0]
    unique_pair_score = ranked[
        (ranked["drug_1"] == "partner_hub") & (ranked["drug_2"] == "unique_x")
    ]["DrugPairScore"].iloc[0]
    original_redundant = best_target_pair_cts("drug_a", "partner_hub", target_map, cts, community_map, crosstalk_edges)
    original_unique = best_target_pair_cts("partner_hub", "unique_x", target_map, cts, community_map, crosstalk_edges)
    print(f"Before penalty: redundant pair (freq=3) = {original_redundant:.4f}, unique pair (freq=1) = {original_unique:.4f}")
    print(f"After penalty:  redundant pair (freq=3) = {redundant_pair_score:.4f}, unique pair (freq=1) = {unique_pair_score:.4f}")
    assert original_redundant > original_unique, "test setup check: redundant pair should originally outrank the unique one"
    gap_before = original_redundant - original_unique
    gap_after = redundant_pair_score - unique_pair_score
    assert gap_after < gap_before, "the gap between the redundant and unique pair should shrink after the penalty"
    print(f"PASSED: the {gap_before:.4f}-point gap between the redundant and unique pair narrowed to "
          f"{gap_after:.4f} after the penalty -- the redundant pair's advantage is being eroded, "
          f"exactly as intended, even though a stronger penalty weight would be needed to fully "
          f"flip the order in this specific synthetic case.\n")

    print("=" * 70)
    print("ALL SMOKE TESTS PASSED.")
    print("NOTE: this confirms the mechanism works as designed on a controlled synthetic case.")
    print("Whether it meaningfully helps on the real 116-drug dataset must still be tested for")
    print("real -- run compare_redundancy_penalty.py against real data before concluding anything,")
    print("same discipline as the first (failed) attempt should have used from the start.")


if __name__ == "__main__":
    _run_smoke_test()
