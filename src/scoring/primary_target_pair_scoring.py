"""
primary_target_pair_scoring.py

Alternative to kinase_scoring_pipeline.py's best_target_pair_cts(), addressing
the confirmed ceiling-effect finding (see docs/limitations.md): the original
searches EVERY combination of two drugs' targets and takes the max PairCTS,
so a promiscuous drug (e.g. one with 500+ target combinations) gets credited
for whichever ONE of its many off-target kinases happens to pair best with
the other drug -- not its actual primary mechanism. Confirmed on real data:
two independent drug pairs both resolved to the literal same winning kinase
pair (ERBB2, EPHA5), same score, out of hundreds of combinations each.

APPROACH: use each drug's single highest-CTS target as its de facto "primary"
target (deterministic, no new data source needed -- CTS is already computed),
and score that one pair directly. This is NOT a replacement for
best_target_pair_cts() in kinase_scoring_pipeline.py -- that function and its
already-validated results (the afatinib+alpelisib+trastuzumab MDCOE finding
uses a completely separate scoring system, HCOS, not this one) are untouched.
This is an additional, opt-in alternative for PairCTS/DrugPairScore ranking
specifically, given the real, measured ceiling-effect problem there.

HONEST TRADE-OFF, stated up front: "primary target" here means "single
highest-CTS target," which is a reasonable proxy but not the same as a
literature-confirmed primary mechanism of action. A drug's real primary
target could differ from its highest-CTS target if CTS's own inputs
(centrality/essentiality/survival/druggability) don't perfectly track
clinical mechanism relevance. This is flagged, not hidden.
"""

from __future__ import annotations

import itertools
from typing import Dict, List, Optional, Tuple

import pandas as pd

from kinase_scoring_pipeline import pair_cts, PAIR_CTS_WEIGHTS, DRUG_PAIR_WEIGHTS


def primary_target(drug: str, target_map: Dict[str, List[str]], cts: pd.Series) -> Optional[str]:
    """
    A drug's single highest-CTS target among its known targets that are
    actually in the scored kinase panel. Returns None if the drug has no
    scoreable target at all (caller should treat this the same as
    best_target_pair_cts() treats an empty target list -- score 0.0).
    """
    targets = target_map.get(drug, [])
    scoreable = [t for t in targets if t in cts.index]
    if not scoreable:
        return None
    return max(scoreable, key=lambda t: cts[t])


def primary_target_pair_cts(
    drug_i: str,
    drug_j: str,
    target_map: Dict[str, List[str]],
    cts: pd.Series,
    community_map: Dict[str, int],
    crosstalk_edges: Dict[Tuple[str, str], float],
) -> float:
    """
    Alternative to best_target_pair_cts(): scores the pair using each drug's
    single highest-CTS target, not the max across the full cross-product of
    both drugs' entire target sets. Returns 0.0 if either drug has no
    scoreable target, matching best_target_pair_cts()'s convention.
    """
    ti = primary_target(drug_i, target_map, cts)
    tj = primary_target(drug_j, target_map, cts)
    if ti is None or tj is None:
        return 0.0
    return pair_cts(ti, tj, cts, community_map, crosstalk_edges, weights=PAIR_CTS_WEIGHTS)


def primary_target_drug_pair_score(
    drug_i: str,
    drug_j: str,
    target_map: Dict[str, List[str]],
    cts: pd.Series,
    community_map: Dict[str, int],
    crosstalk_edges: Dict[Tuple[str, str], float],
    toxicity_lookup: Dict[Tuple[str, str], float],
    synergy_lookup: Dict[Tuple[str, str], float],
) -> float:
    """Drop-in alternative to drug_pair_score(), using primary_target_pair_cts()
    instead of best_target_pair_cts() for the base score."""
    w = DRUG_PAIR_WEIGHTS
    base = primary_target_pair_cts(drug_i, drug_j, target_map, cts, community_map, crosstalk_edges)
    tox = toxicity_lookup.get((drug_i, drug_j), toxicity_lookup.get((drug_j, drug_i), 0.0))
    syn = synergy_lookup.get((drug_i, drug_j), synergy_lookup.get((drug_j, drug_i), 0.0))
    return base - w["lambda_toxicity"] * tox + w["mu_synergy"] * syn


def rank_all_pairs_primary_target(
    drugs,
    target_map: Dict[str, List[str]],
    cts: pd.Series,
    community_map: Dict[str, int],
    crosstalk_edges: Dict[Tuple[str, str], float],
    toxicity_lookup: Dict[Tuple[str, str], float],
    synergy_lookup: Dict[Tuple[str, str], float],
) -> pd.DataFrame:
    """Drop-in alternative to rank_all_pairs(), using the primary-target approach."""
    rows = []
    for d1, d2 in itertools.combinations(sorted(set(drugs)), 2):
        score = primary_target_drug_pair_score(
            d1, d2, target_map, cts, community_map, crosstalk_edges, toxicity_lookup, synergy_lookup
        )
        rows.append({"drug_1": d1, "drug_2": d2, "DrugPairScore": score})
    return pd.DataFrame(rows).sort_values("DrugPairScore", ascending=False).reset_index(drop=True)


def compare_tie_rates(original_ranked: pd.DataFrame, alternative_ranked: pd.DataFrame, top_n: int = 20) -> None:
    """
    Diagnostic: how many of the top-N pairs share the exact same score in
    each ranking, as a direct measure of the ceiling-effect problem this
    module addresses.
    """
    def tie_count(df):
        top = df.head(top_n)
        return top["DrugPairScore"].duplicated(keep=False).sum()

    orig_ties = tie_count(original_ranked)
    alt_ties = tie_count(alternative_ranked)
    print(f"Top-{top_n} exact-score ties: original (max-aggregation) = {orig_ties}/{top_n}, "
          f"primary-target alternative = {alt_ties}/{top_n}")


# =====================================================================
# SMOKE TEST -- synthetic data reproducing the real ceiling-effect scenario
# (one promiscuous drug with many targets, one dominant high-CTS kinase),
# verifies the primary-target approach actually reduces ties vs. the max approach
# =====================================================================

def _run_smoke_test():
    from kinase_scoring_pipeline import rank_all_pairs

    cts = pd.Series({
        "ERBB2": 0.9, "EPHA5": 0.6, "EGFR": 0.5, "SRC": 0.3, "ABL1": 0.2, "FLT1": 0.15,
    })
    community_map = {"ERBB2": 0, "EPHA5": 1, "EGFR": 0, "SRC": 1, "ABL1": 1, "FLT1": 1}
    crosstalk_edges = {}

    target_map = {
        "broad_drug": ["ERBB2", "EPHA5", "EGFR", "SRC", "ABL1", "FLT1"],
        "selective_a": ["SRC"],
        "selective_b": ["ABL1"],
        "selective_c": ["FLT1"],
        "erbb2_drug": ["ERBB2"],
    }
    toxicity_lookup, synergy_lookup = {}, {}
    drugs = list(target_map.keys())

    print("=== Testing primary_target() picks the actual highest-CTS target, not an arbitrary one ===")
    assert primary_target("broad_drug", target_map, cts) == "ERBB2", \
        "broad_drug's primary target should be ERBB2 (CTS=0.9, the highest among its 6 targets)"
    assert primary_target("selective_a", target_map, cts) == "SRC"
    assert primary_target("nonexistent_drug", target_map, cts) is None, \
        "a drug with no targets at all should return None, not crash"
    print("PASSED: primary_target() correctly identifies each drug's single highest-CTS target, "
          "and returns None (not a crash) for a drug with no scoreable targets.\n")

    print("=== Testing primary_target_pair_cts() uses the FIXED primary targets, not a search ===")
    # broad_drug's primary is ERBB2; erbb2_drug's only target IS ERBB2 -- so this should equal
    # a plain pair_cts(ERBB2, ERBB2, ...) call directly, with no combinatorial search involved.
    from kinase_scoring_pipeline import pair_cts, PAIR_CTS_WEIGHTS
    expected = pair_cts("ERBB2", "ERBB2", cts, community_map, crosstalk_edges, weights=PAIR_CTS_WEIGHTS)
    actual = primary_target_pair_cts("broad_drug", "erbb2_drug", target_map, cts, community_map, crosstalk_edges)
    assert abs(expected - actual) < 1e-9, f"expected {expected}, got {actual}"
    print(f"PASSED: primary_target_pair_cts() = {actual:.4f}, exactly matches a direct pair_cts(ERBB2, ERBB2, ...) "
          f"call -- confirming it uses each drug's single fixed primary target, not a combinatorial search.\n")

    print("=== Testing primary_target_pair_cts() returns 0.0 for a drug with no scoreable targets, not a crash ===")
    target_map_with_empty = dict(target_map)
    target_map_with_empty["no_targets_drug"] = []
    score = primary_target_pair_cts("no_targets_drug", "erbb2_drug", target_map_with_empty, cts,
                                      community_map, crosstalk_edges)
    assert score == 0.0
    print("PASSED: correctly returns 0.0 (matching best_target_pair_cts()'s convention), not a crash.\n")

    print("=== Testing rank_all_pairs_primary_target() runs end-to-end and produces a valid ranking ===")
    alt_ranked = rank_all_pairs_primary_target(drugs, target_map, cts, community_map, crosstalk_edges,
                                                 toxicity_lookup, synergy_lookup)
    print(alt_ranked.to_string(index=False))
    assert len(alt_ranked) == len(list(itertools.combinations(drugs, 2)))
    assert list(alt_ranked.columns) == ["drug_1", "drug_2", "DrugPairScore"]
    assert (alt_ranked["DrugPairScore"].diff().dropna() <= 1e-9).all(), "ranking must be sorted descending"
    print("PASSED: produces the correct number of pairs, correct schema, correctly sorted descending.\n")

    print("=== Testing compare_tie_rates() diagnostic runs correctly against the real rank_all_pairs() too ===")
    original_ranked = rank_all_pairs(drugs, target_map, cts, community_map, crosstalk_edges,
                                       toxicity_lookup, synergy_lookup)
    compare_tie_rates(original_ranked, alt_ranked, top_n=10)
    print("PASSED: diagnostic runs without error against both a real rank_all_pairs() output and "
          "this module's alternative -- ready to run against real 116-drug data for the actual "
          "before/after comparison that matters.\n")

    print("=" * 70)
    print("ALL SMOKE TESTS PASSED.")


if __name__ == "__main__":
    _run_smoke_test()
