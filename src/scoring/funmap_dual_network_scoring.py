"""
funmap_dual_network_scoring.py

Implements the three-model ablation experiment specified for evaluating
whether a FunMap-derived redundancy score (Zhang lab, Baylor College of
Medicine; Nature Cancer 2024 -- a pan-cancer functional gene network built
via supervised ML directly on CPTAC proteogenomics data, 10,525 genes,
196,800 associations) contributes real information to PairCTS/TripletCTS
beyond what STRING's crosstalk term already provides.

DESIGN, as specified: a shared weight budget across all three models, so
any difference in results is attributable to the network source, not to
an uncontrolled change in how much total weight the crosstalk term carries.

    Model 1 (baseline):  PairCTS = CTS + STRING           (delta_string=0.10, delta_funmap=0.00)
    Model 2 (FunMap only): PairCTS = CTS + FunMap          (delta_string=0.00, delta_funmap=0.10)
    Model 3 (combined):   PairCTS = CTS + STRING + FunMap  (delta_string=0.05, delta_funmap=0.05)

alpha/beta/gamma (CTS(i), CTS(j), community complementarity) are held fixed
at their original values (0.35/0.35/0.20) in all three models -- only the
crosstalk term's internal composition changes.

WHY NOT ASSUME COMBINED IS BETTER: the real, confirmed CPTAC-vs-STRING
crosstalk cross-check earlier in this project found only 36% overlap
between real measured co-expression and STRING-predicted edges -- meaning
the two networks capture partly different real relationships, not that
either is simply redundant with the other. This experiment tests that
directly rather than assuming it.

NEITHER kinase_scoring_pipeline.py's original pair_cts()/best_target_pair_cts()/
rank_all_pairs() NOR crosstalk_strength()/complementarity() are modified --
this module only adds new functions alongside them.
"""

from __future__ import annotations

import itertools
from typing import Dict, List, Tuple

import pandas as pd

from kinase_scoring_pipeline import complementarity, crosstalk_strength, DRUG_PAIR_WEIGHTS


# =====================================================================
# 1. THE THREE MODELS' WEIGHT SCHEMES (shared budget, as specified)
# =====================================================================

DUAL_NETWORK_MODELS = {
    "model1_string_baseline": {"alpha": 0.35, "beta": 0.35, "gamma": 0.20, "delta_string": 0.10, "delta_funmap": 0.00},
    "model2_funmap_only":     {"alpha": 0.35, "beta": 0.35, "gamma": 0.20, "delta_string": 0.00, "delta_funmap": 0.10},
    "model3_combined":        {"alpha": 0.35, "beta": 0.35, "gamma": 0.20, "delta_string": 0.05, "delta_funmap": 0.05},
}


# =====================================================================
# 2. FUNMAP NETWORK LOADING
# =====================================================================

def load_funmap_network(
    funmap_tsv_path: str,
    kinase_panel: List[str] = None,
) -> Dict[Tuple[str, str], float]:
    """
    Loads a real FunMap network file (gene1, gene2, score -- the format
    documented in the bzhanglab/funmap GitHub repo's output, "funmap.tsv").
    Restricts to the kinase panel if provided (the real FunMap network
    covers 10,525 genes genome-wide; only the ~90 kinases in this project's
    panel are relevant here).

    Gene symbols are uppercased on both sides for matching, since gene
    symbol casing conventions can differ across real data sources (the
    same real issue already found and fixed for drug names earlier in
    this project applies here too -- checked, not assumed, via the real
    file's actual column values before trusting a direct match).
    """
    df = pd.read_csv(funmap_tsv_path, sep="\t")
    expected_cols = {"gene1", "gene2", "score"}
    if not expected_cols.issubset(set(df.columns)):
        raise ValueError(
            f"Expected columns {expected_cols} in {funmap_tsv_path}, got {df.columns.tolist()}. "
            f"Check the real file's header before trusting this loader."
        )

    df["gene1"] = df["gene1"].astype(str).str.upper().str.strip()
    df["gene2"] = df["gene2"].astype(str).str.upper().str.strip()

    if kinase_panel is not None:
        panel_upper = {g.upper() for g in kinase_panel}
        before = len(df)
        df = df[df["gene1"].isin(panel_upper) & df["gene2"].isin(panel_upper)]
        print(f"Restricted FunMap network to kinase panel: {len(df)}/{before} edges retained.")

    edges = {}
    for _, row in df.iterrows():
        edges[(row["gene1"], row["gene2"])] = float(row["score"])
    return edges


# =====================================================================
# 3. DUAL-NETWORK PairCTS / BEST-TARGET-PAIR / RANKING
# =====================================================================

def pair_cts_dual_network(
    kinase_i: str,
    kinase_j: str,
    cts: pd.Series,
    community_map: Dict[str, int],
    string_edges: Dict[Tuple[str, str], float],
    funmap_edges: Dict[Tuple[str, str], float],
    weights: Dict[str, float],
) -> float:
    """
    PairCTS(i,j) = alpha*CTS(i) + beta*CTS(j) + gamma*Complementarity
                   + delta_string*STRING_crosstalk + delta_funmap*FunMap_crosstalk

    Reuses the ORIGINAL, unmodified complementarity() and crosstalk_strength()
    functions from kinase_scoring_pipeline.py for both network lookups --
    crosstalk_strength() is generic over whatever edge dict it's given, so
    no new lookup logic was needed or written for FunMap specifically.
    """
    comp = complementarity(kinase_i, kinase_j, community_map)
    string_cross = crosstalk_strength(kinase_i, kinase_j, string_edges)
    funmap_cross = crosstalk_strength(kinase_i, kinase_j, funmap_edges)
    return (
        weights["alpha"] * cts.get(kinase_i, 0.0)
        + weights["beta"] * cts.get(kinase_j, 0.0)
        + weights["gamma"] * comp
        + weights["delta_string"] * string_cross
        + weights["delta_funmap"] * funmap_cross
    )


def best_target_pair_cts_dual(
    drug_i: str,
    drug_j: str,
    target_map: Dict[str, List[str]],
    cts: pd.Series,
    community_map: Dict[str, int],
    string_edges: Dict[Tuple[str, str], float],
    funmap_edges: Dict[Tuple[str, str], float],
    weights: Dict[str, float],
) -> float:
    """Dual-network analogue of best_target_pair_cts(): same max-over-target-combinations
    search, just scoring each candidate pair with pair_cts_dual_network()."""
    targets_i = target_map.get(drug_i, [])
    targets_j = target_map.get(drug_j, [])
    if not targets_i or not targets_j:
        return 0.0
    return max(
        pair_cts_dual_network(ki, kj, cts, community_map, string_edges, funmap_edges, weights)
        for ki, kj in itertools.product(targets_i, targets_j)
    )


def drug_pair_score_dual(
    drug_i: str,
    drug_j: str,
    target_map: Dict[str, List[str]],
    cts: pd.Series,
    community_map: Dict[str, int],
    string_edges: Dict[Tuple[str, str], float],
    funmap_edges: Dict[Tuple[str, str], float],
    weights: Dict[str, float],
    toxicity_lookup: Dict[Tuple[str, str], float],
    synergy_lookup: Dict[Tuple[str, str], float],
) -> float:
    """Dual-network analogue of drug_pair_score(): identical toxicity/synergy handling,
    identical DRUG_PAIR_WEIGHTS, only the base PairCTS computation differs."""
    w = DRUG_PAIR_WEIGHTS
    base = best_target_pair_cts_dual(drug_i, drug_j, target_map, cts, community_map, string_edges, funmap_edges, weights)
    tox = toxicity_lookup.get((drug_i, drug_j), toxicity_lookup.get((drug_j, drug_i), 0.0))
    syn = synergy_lookup.get((drug_i, drug_j), synergy_lookup.get((drug_j, drug_i), 0.0))
    return base - w["lambda_toxicity"] * tox + w["mu_synergy"] * syn


def rank_all_pairs_dual(
    drugs,
    target_map: Dict[str, List[str]],
    cts: pd.Series,
    community_map: Dict[str, int],
    string_edges: Dict[Tuple[str, str], float],
    funmap_edges: Dict[Tuple[str, str], float],
    weights: Dict[str, float],
    toxicity_lookup: Dict[Tuple[str, str], float],
    synergy_lookup: Dict[Tuple[str, str], float],
) -> pd.DataFrame:
    """Dual-network analogue of rank_all_pairs()."""
    rows = []
    for d1, d2 in itertools.combinations(sorted(set(drugs)), 2):
        score = drug_pair_score_dual(d1, d2, target_map, cts, community_map, string_edges, funmap_edges, weights, toxicity_lookup, synergy_lookup)
        rows.append({"drug_1": d1, "drug_2": d2, "DrugPairScore": score})
    return pd.DataFrame(rows).sort_values("DrugPairScore", ascending=False).reset_index(drop=True)


# =====================================================================
# SMOKE TEST -- confirms Model 1 exactly matches the original pair_cts()
# when delta_funmap=0, confirms Model 2 exactly matches when computed with
# only the FunMap edges swapped in for STRING, and confirms Model 3 is a
# real, correctly-weighted blend of both -- against a synthetic case where
# STRING and FunMap deliberately disagree on one pair, to verify the
# combination logic is doing real arithmetic, not silently favoring one.
# =====================================================================

def _run_smoke_test():
    from kinase_scoring_pipeline import pair_cts, PAIR_CTS_WEIGHTS

    cts = pd.Series({"ERBB2": 0.9, "EPHA5": 0.6, "EGFR": 0.5, "SRC": 0.3})
    community_map = {"ERBB2": 0, "EPHA5": 1, "EGFR": 0, "SRC": 1}
    string_edges = {("ERBB2", "EPHA5"): 0.8}   # STRING says these two are strongly linked
    funmap_edges = {("ERBB2", "EPHA5"): 0.2}   # FunMap disagrees -- weakly linked in real tumors

    print("=== Testing Model 1 (STRING baseline) exactly matches the original pair_cts() ===")
    model1_weights = DUAL_NETWORK_MODELS["model1_string_baseline"]
    dual_score = pair_cts_dual_network("ERBB2", "EPHA5", cts, community_map, string_edges, funmap_edges, model1_weights)
    original_score = pair_cts("ERBB2", "EPHA5", cts, community_map, string_edges, weights=PAIR_CTS_WEIGHTS)
    assert abs(dual_score - original_score) < 1e-12, f"Model 1 should exactly match the original: {dual_score} vs {original_score}"
    print(f"PASSED: Model 1 = {dual_score:.4f}, identical to the original, unmodified pair_cts().\n")

    print("=== Testing Model 2 (FunMap only) uses ONLY the FunMap edge, ignoring STRING entirely ===")
    model2_weights = DUAL_NETWORK_MODELS["model2_funmap_only"]
    model2_score = pair_cts_dual_network("ERBB2", "EPHA5", cts, community_map, string_edges, funmap_edges, model2_weights)
    # Manually compute expected: alpha*0.9 + beta*0.6 + gamma*1.0 (cross-community) + delta_funmap*0.2
    expected_model2 = 0.35*0.9 + 0.35*0.6 + 0.20*1.0 + 0.10*0.2
    assert abs(model2_score - expected_model2) < 1e-9, f"expected {expected_model2}, got {model2_score}"
    # Confirm it does NOT match Model 1 (since STRING vs FunMap edges deliberately disagree)
    assert abs(model2_score - dual_score) > 0.01, "Model 2 should differ from Model 1 given disagreeing edge weights"
    print(f"PASSED: Model 2 = {model2_score:.4f}, correctly uses only the FunMap edge (0.2), "
          f"differs from Model 1 ({dual_score:.4f}) exactly because STRING and FunMap disagree here.\n")

    print("=== Testing Model 3 (combined) is a genuine, correctly-weighted blend of both networks ===")
    model3_weights = DUAL_NETWORK_MODELS["model3_combined"]
    model3_score = pair_cts_dual_network("ERBB2", "EPHA5", cts, community_map, string_edges, funmap_edges, model3_weights)
    expected_model3 = 0.35*0.9 + 0.35*0.6 + 0.20*1.0 + 0.05*0.8 + 0.05*0.2
    assert abs(model3_score - expected_model3) < 1e-9, f"expected {expected_model3}, got {model3_score}"
    # Model 3 should sit strictly between Model 1 (all-STRING) and Model 2 (all-FunMap) on this pair,
    # since it's a genuine blend, not a copy of either
    lo, hi = sorted([dual_score, model2_score])
    assert lo <= model3_score <= hi, f"Model 3 ({model3_score}) should sit between Model 1 and Model 2 when they disagree"
    print(f"PASSED: Model 3 = {model3_score:.4f}, correctly sits between Model 1 ({dual_score:.4f}) and "
          f"Model 2 ({model2_score:.4f}) -- confirming it's a genuine, correctly-weighted blend of both "
          f"networks' real (disagreeing) evidence, not silently defaulting to one.\n")

    print("=== Testing all three weight schemes share the same total budget (fair comparison) ===")
    for name, w in DUAL_NETWORK_MODELS.items():
        total = w["alpha"] + w["beta"] + w["gamma"] + w["delta_string"] + w["delta_funmap"]
        assert abs(total - 1.0) < 1e-9, f"{name} weights sum to {total}, not 1.0 -- comparison would not be fair"
    print("PASSED: all three models' weights sum to exactly 1.0 -- any difference in results is attributable "
          "to the network source, not an uncontrolled change in total weight budget.\n")

    print("=== Testing load_funmap_network() parses a real-format file and restricts to a panel correctly ===")
    import tempfile, os
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write("gene1\tgene2\tscore\n")
        f.write("ERBB2\tEPHA5\t0.734\n")
        f.write("erbb2\tsrc\t0.412\n")  # deliberately lowercase, to test case-normalization
        f.write("BRCA1\tTP53\t0.891\n")  # not in kinase panel -- should be filtered out
        temp_path = f.name
    edges = load_funmap_network(temp_path, kinase_panel=["ERBB2", "EPHA5", "EGFR", "SRC"])
    os.unlink(temp_path)
    assert edges[("ERBB2", "EPHA5")] == 0.734
    assert edges[("ERBB2", "SRC")] == 0.412, "lowercase input should still normalize and match"
    assert ("BRCA1", "TP53") not in edges, "genes outside the kinase panel should be filtered out"
    assert len(edges) == 2
    print("PASSED: loader correctly parses real-format FunMap TSVs, normalizes case, and restricts to "
          "the given kinase panel.\n")

    print("=" * 70)
    print("ALL SMOKE TESTS PASSED.")
    print("Ready to run against real STRING + real FunMap data for the actual three-model comparison.")


if __name__ == "__main__":
    _run_smoke_test()
