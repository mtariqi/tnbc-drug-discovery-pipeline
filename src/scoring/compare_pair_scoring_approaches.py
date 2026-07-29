"""
compare_pair_scoring_approaches.py

Runs BOTH the original max-aggregation PairCTS/DrugPairScore (already in
kinase_scoring_pipeline.py, unmodified) and the new primary-target
alternative (primary_target_pair_scoring.py) against your real data, side by
side, so you can see directly whether the alternative actually reduces the
confirmed ceiling-effect tying and produces a more differentiated ranking --
not just trust that it should in theory.

Run from src/scoring/ (needs run_pairs_and_triplets.py with the drug-name
patch already applied there, kinase_scoring_pipeline.py, and
primary_target_pair_scoring.py all in the same directory):
    python3 compare_pair_scoring_approaches.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from kinase_scoring_pipeline import rank_all_pairs
from primary_target_pair_scoring import rank_all_pairs_primary_target, compare_tie_rates

try:
    from run_pairs_and_triplets import load_real_inputs
except ImportError:
    raise SystemExit(
        "run_pairs_and_triplets.py (with the drug-name patch already applied) must be "
        "in the same directory -- copy it here first."
    )


def build_toxicity_lookup_from_local_faers(faers_dir, drug_list):
    """Same logic as run_pairs_and_triplets.py's version, duplicated here to avoid an
    import-order dependency -- keeps this script runnable on its own."""
    import pandas as pd
    faers_path = Path(faers_dir)
    reaction_sets = {}
    for drug in drug_list:
        f = faers_path / f"{drug}_faers.tsv"
        if f.exists():
            df = pd.read_csv(f, sep="\t")
            reaction_sets[drug] = set(df["reaction"].str.lower()) if "reaction" in df.columns else set()
    toxicity_lookup = {}
    import itertools
    for d1, d2 in itertools.combinations(drug_list, 2):
        if d1 in reaction_sets and d2 in reaction_sets:
            shared = reaction_sets[d1] & reaction_sets[d2]
            total = reaction_sets[d1] | reaction_sets[d2]
            toxicity_lookup[(d1, d2)] = len(shared) / len(total) if total else 0.0
    return toxicity_lookup


def main():
    base_dir = sys.argv[1] if len(sys.argv) > 1 else str(Path.home() / "rtk_nrtk_tnbc")
    drug_list_path = sys.argv[2] if len(sys.argv) > 2 else str(Path(base_dir) / "data/raw/drugs/drug_list.txt")

    from run_pairs_and_triplets import _normalize_drug_name
    drug_list = [_normalize_drug_name(d.strip()) for d in open(drug_list_path) if d.strip()]

    print("Loading real CTS, community, crosstalk, and DGIdb data...")
    cts, community_map, crosstalk_edges, target_map = load_real_inputs(base_dir)

    scoreable_drugs = [d for d in drug_list if d in target_map and any(k in cts.index for k in target_map[d])]
    print(f"  {len(scoreable_drugs)}/{len(drug_list)} drugs scoreable")

    # Same 60-drug cap as run_pairs_and_triplets.py, for a fair, comparably-sized comparison
    max_drugs = 60
    if len(scoreable_drugs) > max_drugs:
        drug_max_cts = {d: max(cts.get(k, 0) for k in target_map[d]) for d in scoreable_drugs}
        scoreable_drugs = sorted(scoreable_drugs, key=lambda d: drug_max_cts[d], reverse=True)[:max_drugs]
        print(f"  Capped to top {max_drugs} drugs by max target CTS")

    toxicity_lookup = build_toxicity_lookup_from_local_faers(
        str(Path(base_dir) / "data/processed/faers"), scoreable_drugs
    )
    synergy_lookup = {}

    print("\nRunning ORIGINAL (max-aggregation) ranking...")
    original_ranked = rank_all_pairs(scoreable_drugs, target_map, cts, community_map, crosstalk_edges,
                                       toxicity_lookup, synergy_lookup)

    print("Running PRIMARY-TARGET alternative ranking...")
    alt_ranked = rank_all_pairs_primary_target(scoreable_drugs, target_map, cts, community_map, crosstalk_edges,
                                                 toxicity_lookup, synergy_lookup)

    print("\n=== Top 20: ORIGINAL (max-aggregation) ===")
    print(original_ranked.head(20).to_string(index=False))

    print("\n=== Top 20: PRIMARY-TARGET alternative ===")
    print(alt_ranked.head(20).to_string(index=False))

    print()
    compare_tie_rates(original_ranked, alt_ranked, top_n=20)

    n_unique_scores_orig = original_ranked["DrugPairScore"].nunique()
    n_unique_scores_alt = alt_ranked["DrugPairScore"].nunique()
    print(f"\nUnique scores across the FULL ranking ({len(original_ranked)} pairs): "
          f"original = {n_unique_scores_orig}, primary-target = {n_unique_scores_alt}")

    out_path = Path(base_dir) / "data/processed/pair_ranking_primary_target_alternative.tsv"
    alt_ranked.to_csv(out_path, sep="\t", index=False)
    print(f"\nWrote primary-target ranking to {out_path} for your own further review.")


if __name__ == "__main__":
    main()
