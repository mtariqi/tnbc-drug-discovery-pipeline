"""
compare_redundancy_penalty.py

Runs the original max-aggregation PairCTS ranking and the new redundancy-
penalized alternative against your real data, side by side. Same discipline
as compare_pair_scoring_approaches.py (used for the first, rejected attempt):
this measures the real effect, it doesn't assume one.

Run from src/scoring/:
    python3 compare_redundancy_penalty.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from kinase_scoring_pipeline import rank_all_pairs
from redundancy_penalized_pair_scoring import rank_all_pairs_redundancy_penalized

try:
    from run_pairs_and_triplets import load_real_inputs, _normalize_drug_name
except ImportError:
    raise SystemExit(
        "run_pairs_and_triplets.py (with the drug-name patch already applied) must be "
        "in the same directory -- copy it here first."
    )


def build_toxicity_lookup_from_local_faers(faers_dir, drug_list):
    import pandas as pd
    import itertools
    faers_path = Path(faers_dir)
    reaction_sets = {}
    for drug in drug_list:
        f = faers_path / f"{drug}_faers.tsv"
        if f.exists():
            df = pd.read_csv(f, sep="\t")
            reaction_sets[drug] = set(df["reaction"].str.lower()) if "reaction" in df.columns else set()
    toxicity_lookup = {}
    for d1, d2 in itertools.combinations(drug_list, 2):
        if d1 in reaction_sets and d2 in reaction_sets:
            shared = reaction_sets[d1] & reaction_sets[d2]
            total = reaction_sets[d1] | reaction_sets[d2]
            toxicity_lookup[(d1, d2)] = len(shared) / len(total) if total else 0.0
    return toxicity_lookup


def main():
    base_dir = sys.argv[1] if len(sys.argv) > 1 else str(Path.home() / "rtk_nrtk_tnbc")
    penalty_weight = float(sys.argv[2]) if len(sys.argv) > 2 else 0.05
    drug_list_path = str(Path(base_dir) / "data/raw/drugs/drug_list.txt")

    drug_list = [_normalize_drug_name(d.strip()) for d in open(drug_list_path) if d.strip()]

    print("Loading real CTS, community, crosstalk, and DGIdb data...")
    cts, community_map, crosstalk_edges, target_map = load_real_inputs(base_dir)

    scoreable_drugs = [d for d in drug_list if d in target_map and any(k in cts.index for k in target_map[d])]
    print(f"  {len(scoreable_drugs)}/{len(drug_list)} drugs scoreable")

    max_drugs = 60
    if len(scoreable_drugs) > max_drugs:
        drug_max_cts = {d: max(cts.get(k, 0) for k in target_map[d]) for d in scoreable_drugs}
        scoreable_drugs = sorted(scoreable_drugs, key=lambda d: drug_max_cts[d], reverse=True)[:max_drugs]
        print(f"  Capped to top {max_drugs} drugs by max target CTS")

    toxicity_lookup = build_toxicity_lookup_from_local_faers(
        str(Path(base_dir) / "data/processed/faers"), scoreable_drugs
    )
    synergy_lookup = {}

    print(f"\nRunning ORIGINAL (max-aggregation) ranking...")
    original_ranked = rank_all_pairs(scoreable_drugs, target_map, cts, community_map, crosstalk_edges,
                                       toxicity_lookup, synergy_lookup)

    print(f"Running REDUNDANCY-PENALIZED ranking (penalty_weight={penalty_weight})...")
    penalized_ranked = rank_all_pairs_redundancy_penalized(
        scoreable_drugs, target_map, cts, community_map, crosstalk_edges,
        toxicity_lookup, synergy_lookup, redundancy_penalty_weight=penalty_weight,
    )

    print("\n=== Top 20: ORIGINAL (max-aggregation) ===")
    print(original_ranked.head(20).to_string(index=False))

    print("\n=== Top 20: REDUNDANCY-PENALIZED ===")
    print(penalized_ranked.head(20)[["drug_1", "drug_2", "DrugPairScore", "winning_kinase_pair", "winning_pair_frequency"]].to_string(index=False))

    n_unique_original = original_ranked["DrugPairScore"].nunique()
    n_unique_penalized = penalized_ranked["DrugPairScore"].nunique()
    print(f"\nUnique scores across the FULL ranking ({len(original_ranked)} pairs): "
          f"original = {n_unique_original}, redundancy-penalized = {n_unique_penalized}")

    top20_ties_original = original_ranked.head(20)["DrugPairScore"].duplicated(keep=False).sum()
    top20_ties_penalized = penalized_ranked.head(20)["DrugPairScore"].duplicated(keep=False).sum()
    print(f"Top-20 exact-score ties: original = {top20_ties_original}/20, "
          f"redundancy-penalized = {top20_ties_penalized}/20")

    # How many DIFFERENT winning kinase pairs appear in the top 20, before vs after --
    # a more direct measure of "is the top of the ranking more diverse now"
    print(f"\nDistinct winning kinase pairs in top 20 (redundancy-penalized): "
          f"{penalized_ranked.head(20)['winning_kinase_pair'].nunique()}")

    out_path = Path(base_dir) / "data/processed/pair_ranking_redundancy_penalized.tsv"
    penalized_ranked.to_csv(out_path, sep="\t", index=False)
    print(f"\nWrote redundancy-penalized ranking to {out_path} for your own further review.")

    print("\nHONEST NOTE: this compares ONE penalty_weight value (0.05 by default). If the effect "
          "looks weak, try re-running with a larger weight, e.g.:")
    print("  python3 compare_redundancy_penalty.py ~/rtk_nrtk_tnbc 0.15")


if __name__ == "__main__":
    main()
