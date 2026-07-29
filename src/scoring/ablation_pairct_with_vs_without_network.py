"""
Ablation test: does adding TNBC-correlation-network-derived crosstalk/
compensation information change PairCTS rankings in a way that matters?

This directly implements the falsifiable test proposed for evaluating whether
network-based redundancy information is worth adding to CTS: rank all drug
pairs twice (STRING-only crosstalk vs. STRING + correlation-network crosstalk
combined), then check whether the manuscript's actual regimens of interest
(trastuzumab+capivasertib, trastuzumab+alpelisib, afatinib-based combinations)
move meaningfully. If they do not move, per the pre-registered criterion
already agreed on, this network layer should be deprioritized rather than
added for its own sake.

This script does NOT modify pair_cts(), best_target_pair_cts(), or
rank_all_pairs() -- it calls them exactly as they already exist, twice, with
two different crosstalk_edges dicts. The correlation-network's crosstalk_edges
(from tnbc_correlation_network.py) uses |r| as edge weight, combined
(dict-merged) with the existing STRING-derived edges for the "with network"
condition; STRING alone for the "without" baseline condition.

Usage (once you have the real pipeline artifacts -- see REQUIRED INPUTS below):
  python -m src.scoring.ablation_pairct_with_vs_without_network \
      --string-crosstalk path/to/string_crosstalk.json \
      --correlation-network results/tables/tnbc_correlation_network.json \
      --cts path/to/cts_scores.csv \
      --community-map path/to/string_communities.json \
      --target-map path/to/drug_target_map.json \
      --toxicity path/to/toxicity_lookup.json \
      --synergy path/to/synergy_lookup.json \
      --focal-pairs "trastuzumab,capivasertib" "trastuzumab,alpelisib" \
      --output results/tables/ablation_network_vs_no_network.csv

REQUIRED INPUTS (all should already exist from the real CTS/PairCTS pipeline
-- this script does not recompute CTS, target maps, or STRING communities
itself, since those are already built and validated elsewhere in this project):
  --string-crosstalk : JSON {"geneA|geneB": weight, ...} of existing STRING edges
  --cts              : CSV with columns gene,cts_score
  --community-map    : JSON {"gene": community_id, ...} from STRING Louvain
  --target-map       : JSON {"drug_name": ["gene1", "gene2", ...], ...}
  --toxicity         : JSON {"drugA|drugB": toxicity_value, ...}
  --synergy          : JSON {"drugA|drugB": synergy_value, ...}
"""
from __future__ import annotations

import argparse
import json

import pandas as pd

from src.scoring.kinase_scoring_pipeline import rank_all_pairs


def _load_edge_json(path: str) -> dict:
    with open(path) as f:
        raw = json.load(f)
    return {tuple(k.split("|")): v for k, v in raw.items()}


def _find_pair_rank(ranked_df: pd.DataFrame, drug_a: str, drug_b: str) -> dict:
    """Rank is 1-indexed position in the sorted DataFrame; order-insensitive
    (checks both drug_1/drug_2 orderings, since rank_all_pairs sorts drug
    names alphabetically per pair, not by input order)."""
    match = ranked_df[
        ((ranked_df["drug_1"] == drug_a) & (ranked_df["drug_2"] == drug_b))
        | ((ranked_df["drug_1"] == drug_b) & (ranked_df["drug_2"] == drug_a))
    ]
    if match.empty:
        return {"found": False, "rank": None, "score": None}
    idx = match.index[0]
    return {"found": True, "rank": int(idx) + 1, "score": float(match.iloc[0]["DrugPairScore"])}


def run_ablation(
    drugs: list[str],
    target_map: dict,
    cts: pd.Series,
    community_map: dict,
    string_crosstalk: dict,
    correlation_crosstalk: dict,
    toxicity_lookup: dict,
    synergy_lookup: dict,
    focal_pairs: list[tuple[str, str]],
) -> pd.DataFrame:
    combined_crosstalk = {**string_crosstalk, **correlation_crosstalk}

    print("Ranking all pairs -- STRING crosstalk only (baseline)...")
    ranking_without = rank_all_pairs(drugs, target_map, cts, community_map, string_crosstalk, toxicity_lookup, synergy_lookup)

    print("Ranking all pairs -- STRING + TNBC correlation-network crosstalk combined...")
    ranking_with = rank_all_pairs(drugs, target_map, cts, community_map, combined_crosstalk, toxicity_lookup, synergy_lookup)

    rows = []
    for drug_a, drug_b in focal_pairs:
        without = _find_pair_rank(ranking_without, drug_a, drug_b)
        with_ = _find_pair_rank(ranking_with, drug_a, drug_b)
        rank_change = (without["rank"] - with_["rank"]) if (without["found"] and with_["found"]) else None
        rows.append({
            "drug_pair": f"{drug_a}+{drug_b}",
            "rank_without_network": without["rank"],
            "rank_with_network": with_["rank"],
            "rank_change": rank_change,  # positive = moved UP (better) with network info
            "score_without_network": without["score"],
            "score_with_network": with_["score"],
        })

    result = pd.DataFrame(rows)
    n_total_pairs = len(ranking_without)
    print(f"\n{n_total_pairs} total pairs ranked. Focal-pair results:")
    print(result.to_string(index=False))
    print(
        "\nInterpretation guide (per the pre-registered criterion): if focal pairs move "
        "UP meaningfully (rank_change clearly positive, not just +/-1-2 out of "
        f"{n_total_pairs} pairs) because they connect highly-connected or compensatory "
        "network regions, the network layer adds real value. If rank_change is small or "
        "near zero across the board, deprioritize this component rather than add it for "
        "its own sake."
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--string-crosstalk", required=True)
    parser.add_argument("--correlation-network", required=True,
                         help="Output of tnbc_correlation_network.py")
    parser.add_argument("--cts", required=True)
    parser.add_argument("--community-map", required=True)
    parser.add_argument("--target-map", required=True)
    parser.add_argument("--toxicity", required=True)
    parser.add_argument("--synergy", required=True)
    parser.add_argument("--focal-pairs", nargs="+", required=True,
                         help='e.g. "trastuzumab,capivasertib" "trastuzumab,alpelisib"')
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    string_crosstalk = _load_edge_json(args.string_crosstalk)

    with open(args.correlation_network) as f:
        corr_net = json.load(f)
    correlation_crosstalk = {tuple(k.split("|")): v for k, v in corr_net["crosstalk_edges"].items()}

    cts = pd.read_csv(args.cts, index_col=0).iloc[:, 0]
    with open(args.community_map) as f:
        community_map = json.load(f)
    with open(args.target_map) as f:
        target_map = json.load(f)
    toxicity_lookup = _load_edge_json(args.toxicity)
    synergy_lookup = _load_edge_json(args.synergy)

    drugs = list(target_map.keys())
    focal_pairs = [tuple(p.split(",")) for p in args.focal_pairs]

    result = run_ablation(
        drugs, target_map, cts, community_map, string_crosstalk, correlation_crosstalk,
        toxicity_lookup, synergy_lookup, focal_pairs,
    )
    result.to_csv(args.output, index=False)
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
