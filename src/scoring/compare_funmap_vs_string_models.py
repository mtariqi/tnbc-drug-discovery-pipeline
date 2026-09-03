"""
compare_funmap_vs_string_models.py

Runs the three-model ablation experiment exactly as specified:
    Model 1 (baseline):   PairCTS = CTS + STRING
    Model 2 (FunMap only): PairCTS = CTS + FunMap           (STRING removed)
    Model 3 (combined):   PairCTS = CTS + STRING + FunMap   (combined)

All three share the same total weight budget (alpha/beta/gamma fixed,
delta_string + delta_funmap always summing to the original delta=0.10),
so any difference in results is attributable to the network source, not
an uncontrolled change in how much weight the crosstalk term carries.

Compares, as specified:
    - top-ranked pairs (overlap between models)
    - top-ranked triplets -- NOT run here; see the module docstring in
      funmap_dual_network_scoring.py for why this isn't currently
      well-defined (TripletCTS has no network-crosstalk term to swap).
    - rank stability (Spearman correlation between full rankings)
    - PairCTS score distribution (checking for score inflation, not just
      score change -- a real improvement should look like re-ranking
      based on new information, not every score simply going up)

HONEST FRAMING: this does not assume Model 3 is better. The real,
confirmed CPTAC-vs-STRING crosstalk cross-check earlier in this project
found only 36% overlap between real measured co-expression and
STRING-predicted edges -- meaning the two networks likely capture partly
different real relationships. This script measures that directly.

Run from src/scoring/ (needs run_pairs_and_triplets.py, with the drug-name
patch already applied, and funmap_dual_network_scoring.py in the same
directory):
    python3 compare_funmap_vs_string_models.py <path_to_real_funmap.tsv>
"""

import sys
from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).parent))

from funmap_dual_network_scoring import DUAL_NETWORK_MODELS, load_funmap_network, rank_all_pairs_dual
from redundancy_penalized_dual_network_scoring import rank_all_pairs_dual_redundancy_penalized

try:
    from run_pairs_and_triplets import load_real_inputs, _normalize_drug_name
except ImportError:
    raise SystemExit(
        "run_pairs_and_triplets.py (with the drug-name patch already applied) must be "
        "in the same directory -- copy it here first."
    )


def build_toxicity_lookup_from_local_faers(faers_dir, drug_list):
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
    if len(sys.argv) < 2:
        raise SystemExit(
            "Usage: python3 compare_funmap_vs_string_models.py <path_to_real_funmap_network.tsv> [base_dir]\n"
            "The FunMap TSV must have columns: gene1, gene2, score (the real bzhanglab/funmap output format)."
        )
    funmap_path = sys.argv[1]
    base_dir = sys.argv[2] if len(sys.argv) > 2 else str(Path.home() / "rtk_nrtk_tnbc")

    print("Loading real CTS, community, STRING crosstalk, and DGIdb data...")
    cts, community_map, string_edges, target_map = load_real_inputs(base_dir)
    print(f"  CTS: {len(cts)} kinases, {len(string_edges)} STRING edges")

    kinase_panel = list(cts.index)
    print(f"\nLoading real FunMap network from {funmap_path}...")
    funmap_edges = load_funmap_network(funmap_path, kinase_panel=kinase_panel)
    print(f"  {len(funmap_edges)} FunMap edges within the {len(kinase_panel)}-kinase panel")

    n_string_pairs = set(tuple(sorted(k)) for k in string_edges.keys())
    n_funmap_pairs = set(tuple(sorted(k)) for k in funmap_edges.keys())
    overlap = n_string_pairs & n_funmap_pairs
    print(f"\n  STRING-FunMap edge overlap (within kinase panel): {len(overlap)}/{len(n_string_pairs | n_funmap_pairs)} "
          f"({100*len(overlap)/max(len(n_string_pairs | n_funmap_pairs),1):.0f}% of the union)")

    drug_list_path = str(Path(base_dir) / "data/raw/drugs/drug_list.txt")
    drug_list = [_normalize_drug_name(d.strip()) for d in open(drug_list_path) if d.strip()]
    scoreable_drugs = [d for d in drug_list if d in target_map and any(k in cts.index for k in target_map[d])]
    max_drugs = 60
    if len(scoreable_drugs) > max_drugs:
        drug_max_cts = {d: max(cts.get(k, 0) for k in target_map[d]) for d in scoreable_drugs}
        scoreable_drugs = sorted(scoreable_drugs, key=lambda d: drug_max_cts[d], reverse=True)[:max_drugs]
    print(f"\n  {len(scoreable_drugs)} drugs scoreable (capped at {max_drugs} by max target CTS)")

    toxicity_lookup = build_toxicity_lookup_from_local_faers(str(Path(base_dir) / "data/processed/faers"), scoreable_drugs)
    synergy_lookup = {}

    faers_dir = Path(base_dir) / "data/processed/faers"
    drugs_with_faers = sum(1 for d in scoreable_drugs if (faers_dir / f"{d}_faers.tsv").exists())
    total_possible_pairs = len(scoreable_drugs) * (len(scoreable_drugs) - 1) // 2
    print(f"\n  FAERS coverage: {drugs_with_faers}/{len(scoreable_drugs)} scoreable drugs have a real local FAERS file "
          f"({100*drugs_with_faers/max(len(scoreable_drugs),1):.0f}%)")
    print(f"  Toxicity-supported pairs: {len(toxicity_lookup)}/{total_possible_pairs} "
          f"({100*len(toxicity_lookup)/max(total_possible_pairs,1):.0f}%) -- pairs without both drugs' FAERS data "
          f"get a toxicity penalty of exactly 0.0, not a missing value; treat this as a real limitation, not neutral.")

    rankings = {}
    rankings_raw = {}
    for model_name, weights in DUAL_NETWORK_MODELS.items():
        print(f"\nRanking pairs under {model_name} (weights={weights})...")
        rankings[model_name] = rank_all_pairs_dual_redundancy_penalized(
            scoreable_drugs, target_map, cts, community_map, string_edges, funmap_edges, weights,
            toxicity_lookup, synergy_lookup, redundancy_penalty_weight=0.15,
        )
        rankings_raw[model_name] = rank_all_pairs_dual(
            scoreable_drugs, target_map, cts, community_map, string_edges, funmap_edges, weights,
            toxicity_lookup, synergy_lookup,
        )
    print("\nNOTE: the redundancy penalty (weight=0.15) is applied ONLY for the top-N tables and top-K overlap")
    print("statistics below, matching how it is validated and adopted in the manuscript (Section 3.8/3.9 --")
    print("'adopted as the default for top-N reporting only'). Real-data testing found this penalty, when applied")
    print("across the FULL ranking rather than restricted to top-N, produces penalties 15-35x larger than the base")
    print("scores themselves for hub-heavy kinases (e.g. ERBB2), pushing over half of all pairs' scores negative --")
    print("an over-penalization artifact, not a real finding. Full-ranking statistics (Spearman, rank-shift, score")
    print("distribution) below therefore use the UNPENALIZED base ranking instead.")

    print("\n" + "=" * 70)
    print("=== Top 15 pairs under each model ===")
    for model_name, df in rankings.items():
        print(f"\n--- {model_name} ---")
        print(df.head(15).to_string(index=False))

    m1, m2, m3 = "model1_string_baseline", "model2_funmap_only", "model3_combined"

    def jaccard_at_k(df_a, df_b, k):
        set_a = set(zip(df_a.head(k)["drug_1"], df_a.head(k)["drug_2"]))
        set_b = set(zip(df_b.head(k)["drug_1"], df_b.head(k)["drug_2"]))
        union = set_a | set_b
        return len(set_a & set_b), len(union), (len(set_a & set_b) / len(union) if union else 0.0)

    print("\n" + "=" * 70)
    print("=== Top-K overlap between models (raw count and Jaccard index, at K=10/20/50) ===")
    pairs_to_compare = [("Model 1 vs Model 2 (STRING vs FunMap only)", m1, m2),
                        ("Model 1 vs Model 3 (baseline vs combined)", m1, m3),
                        ("Model 2 vs Model 3 (FunMap vs combined)", m2, m3)]
    for label, a, b in pairs_to_compare:
        for k in (10, 20, 50):
            n_common, n_union, jacc = jaccard_at_k(rankings[a], rankings[b], k)
            print(f"{label}: top-{k} -> {n_common}/{k} common, Jaccard = {jacc:.3f}")

    print("\n=== Rank stability on the UNPENALIZED base ranking (explicit ranks; see note above on why the") 
    print("    redundancy-penalized ranking is not used for this full-ranking statistic) ===")
    merged = rankings_raw[m1][["drug_1", "drug_2", "DrugPairScore"]].merge(
        rankings_raw[m2][["drug_1", "drug_2", "DrugPairScore"]], on=["drug_1", "drug_2"], suffixes=("_m1", "_m2")
    ).merge(
        rankings_raw[m3][["drug_1", "drug_2", "DrugPairScore"]].rename(columns={"DrugPairScore": "DrugPairScore_m3"}),
        on=["drug_1", "drug_2"],
    )
    merged["rank_m1"] = merged["DrugPairScore_m1"].rank(ascending=False, method="average")
    merged["rank_m2"] = merged["DrugPairScore_m2"].rank(ascending=False, method="average")
    merged["rank_m3"] = merged["DrugPairScore_m3"].rank(ascending=False, method="average")

    rho_12, _ = spearmanr(merged["rank_m1"], merged["rank_m2"])
    rho_13, _ = spearmanr(merged["rank_m1"], merged["rank_m3"])
    rho_23, _ = spearmanr(merged["rank_m2"], merged["rank_m3"])
    print(f"Model 1 vs Model 2: Spearman rho = {rho_12:.3f}")
    print(f"Model 1 vs Model 3: Spearman rho = {rho_13:.3f}")
    print(f"Model 2 vs Model 3: Spearman rho = {rho_23:.3f}")
    print("(rho close to 1.0 = rankings barely changed; lower = the network source is genuinely reordering candidates)")

    print("\n=== Rank-shift magnitude (more directly interpretable than correlation alone) ===")
    for label, a_key, b_key in [("Model 1 vs Model 2", "rank_m1", "rank_m2"),
                                  ("Model 1 vs Model 3", "rank_m1", "rank_m3"),
                                  ("Model 2 vs Model 3", "rank_m2", "rank_m3")]:
        abs_shift = (merged[a_key] - merged[b_key]).abs()
        n_total = len(merged)
        n_moved_50 = (abs_shift > 50).sum()
        print(f"{label}: median absolute rank shift = {abs_shift.median():.1f} positions "
              f"(out of {n_total} ranked pairs); {n_moved_50} pairs ({100*n_moved_50/n_total:.1f}%) moved >50 positions")

    print("\n=== New entrants to top-20 (in B's top-20 but not A's) ===")
    for label, a, b in pairs_to_compare:
        set_a20 = set(zip(rankings[a].head(20)["drug_1"], rankings[a].head(20)["drug_2"]))
        set_b20 = set(zip(rankings[b].head(20)["drug_1"], rankings[b].head(20)["drug_2"]))
        new_in_b = set_b20 - set_a20
        print(f"{label}: {len(new_in_b)} pairs in the second model's top-20 were not in the first's")

    print("\n=== Hub-bias diagnostic: is a small number of kinases dominating the winning pairs? ===")
    for model_name, df in rankings.items():
        kinase_counts = {}
        for wp in df["winning_kinase_pair"].dropna():
            for kinase in wp:
                kinase_counts[kinase] = kinase_counts.get(kinase, 0) + 1
        total_appearances = sum(kinase_counts.values())
        top5 = sorted(kinase_counts.items(), key=lambda x: -x[1])[:5]
        top1_share = (top5[0][1] / total_appearances * 100) if top5 and total_appearances else 0.0
        top5_share = (sum(c for _, c in top5) / total_appearances * 100) if total_appearances else 0.0
        print(f"{model_name}: top kinase '{top5[0][0] if top5 else 'n/a'}' accounts for {top1_share:.1f}% of all "
              f"winning-pair appearances; top-5 kinases account for {top5_share:.1f}%. Top 5: {top5}")
    print("(if one or two kinases dominate every model similarly, that may reflect genuine biological hub status --")
    print(" e.g. ERBB2/EGFR -- rather than a network artifact; compare against each kinase's real STRING/FunMap degree")
    print(" before concluding either way.)")

    print("\n=== PairCTS score distribution per model, UNPENALIZED base ranking (checking for score inflation, not just score change) ===")
    for model_name, df in rankings_raw.items():
        s = df["DrugPairScore"]
        print(f"{model_name}: mean={s.mean():.4f}, std={s.std():.4f}, min={s.min():.4f}, max={s.max():.4f}")
    print("(if Model 3's mean is simply higher than Model 1's with similar spread, that suggests additive score")
    print(" inflation rather than genuine re-ranking -- check the top-20 overlap and rank-stability numbers above")
    print(" together with this, not in isolation.)")

    print("\n" + "=" * 70)
    print("=== Interpretation guide (check the real numbers above against these patterns) ===")
    print("Scenario A -- Top-20 overlap 18-20/20, Spearman >0.95: FunMap provides little new information vs. STRING.")
    print("Scenario B -- Top-20 overlap 8-15/20, Spearman 0.7-0.9: FunMap meaningfully reorders candidates")
    print("              (the most biologically interesting outcome).")
    print("Scenario C -- Combined model's mean score much higher, Spearman very high, top-20 nearly unchanged:")
    print("              score inflation rather than new biological signal -- re-check the distribution stats above.")
    print("Scenario D -- STRING-FunMap edge overlap low, Spearman low, top-20 overlap low: the two networks encode")
    print("              substantially different biology -- consistent with the CPTAC-vs-STRING 36% overlap already")
    print("              found in Section 3.6 of the manuscript.")

    out_dir = Path(base_dir) / "data/processed"
    for model_name, df in rankings.items():
        out_path = out_dir / f"pair_ranking_{model_name}_topN_penalized.tsv"
        df.to_csv(out_path, sep="\t", index=False)
    for model_name, df in rankings_raw.items():
        out_path = out_dir / f"pair_ranking_{model_name}_unpenalized.tsv"
        df.to_csv(out_path, sep="\t", index=False)
    print(f"\nWrote redundancy-penalized rankings (used for top-N/Jaccard stats) to {out_dir}/pair_ranking_model*_topN_penalized.tsv")
    print(f"Wrote unpenalized rankings (used for full-ranking stats) to {out_dir}/pair_ranking_model*_unpenalized.tsv")

    print("\nNOTE: top-ranked TRIPLET comparison is not run here -- TripletCTS has no network-crosstalk term")
    print("in the current codebase to swap STRING for FunMap in. See funmap_dual_network_scoring.py's module")
    print("docstring for what a triplet-level extension would require, and confirm the design before building it.")


if __name__ == "__main__":
    main()
