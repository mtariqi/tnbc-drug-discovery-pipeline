"""
Analysis 1: Weight-perturbation robustness for CTS and PairCTS.

Real inputs only: kinase_scoring_pipeline.compute_cts()/pair_cts() (unmodified),
CTS_all_90_kinases.tsv (real, verified to reproduce published CTS values exactly),
string_edges.tsv (real, 464 edges), community_map.tsv (real, 90 kinases).
"""
from __future__ import annotations

import itertools
import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from kinase_scoring_pipeline import compute_cts, pair_cts, CTS_WEIGHTS, PAIR_CTS_WEIGHTS

RNG_SEED = 42

# Matches this project's standard layout (see README.md "Repository layout").
# Override by setting the TNBC_BASE_DIR environment variable if your checkout
# lives somewhere other than ~/rtk_nrtk_tnbc.
# Resolves to the repo root this script actually lives in (src/scoring/../.. ->
# repo root), so this works regardless of what the repo is named or cloned to,
# rather than guessing a fixed path like ~/rtk_nrtk_tnbc. Override with the
# TNBC_BASE_DIR environment variable if you need to point elsewhere (e.g. running
# this script from outside the repo, or against a different repo's data).
_DEFAULT_BASE_DIR = Path(__file__).resolve().parent.parent.parent
BASE_DIR = Path(os.environ.get("TNBC_BASE_DIR", _DEFAULT_BASE_DIR))
OUTPUT_DIR = BASE_DIR / "results" / "tables" / "robustness_ablation_null"


def load_real_inputs():
    raw = pd.read_csv(BASE_DIR / "data/processed/CTS_all_90_kinases.tsv", sep="\t", index_col="kinase_id")
    kinase_df = raw[["betweenness", "pagerank", "degree", "depmap_score", "cox_hr", "logrank_p",
                      "dgidb_score", "chembl_count", "trial_stage"]].copy()
    edges = pd.read_csv(BASE_DIR / "data/processed/string/string_edges.tsv", sep="\t")
    crosstalk_edges = {(r.source, r.target): r.weight for r in edges.itertuples(index=False)}
    comm = pd.read_csv(BASE_DIR / "data/processed/string/community_map.tsv", sep="\t")
    community_map = dict(zip(comm["kinase_id"], comm["community"]))
    return kinase_df, crosstalk_edges, community_map


def top_k_overlap(rank_a: pd.Series, rank_b: pd.Series, k: int) -> float:
    top_a = set(rank_a.sort_values(ascending=False).head(k).index)
    top_b = set(rank_b.sort_values(ascending=False).head(k).index)
    return len(top_a & top_b) / k


def cts_with_weights(kinase_df, weights):
    return compute_cts(kinase_df, weights=weights)["CTS"]


def pairct_with_weights(cts, community_map, crosstalk_edges, weights, kinases):
    rows = []
    for ki, kj in itertools.combinations(kinases, 2):
        rows.append(((ki, kj), pair_cts(ki, kj, cts, community_map, crosstalk_edges, weights=weights)))
    return pd.Series({pair: score for pair, score in rows})


def normalize_weights(w: dict) -> dict:
    total = sum(w.values())
    return {k: v / total for k, v in w.items()}


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    kinase_df, crosstalk_edges, community_map = load_real_inputs()
    kinases = list(kinase_df.index)

    baseline_cts = cts_with_weights(kinase_df, CTS_WEIGHTS)
    baseline_paircts = pairct_with_weights(baseline_cts, community_map, crosstalk_edges, PAIR_CTS_WEIGHTS, kinases)

    print(f"Baseline CTS computed for {len(baseline_cts)} kinases")
    print(f"Baseline PairCTS computed for {len(baseline_paircts)} pairs")

    # ---------------- Scenario A: +/-25% one-at-a-time perturbations (CTS) ----------------
    scenario_a_results = []
    for component in CTS_WEIGHTS:
        for frac in (-0.25, 0.25):
            perturbed = dict(CTS_WEIGHTS)
            perturbed[component] = CTS_WEIGHTS[component] * (1 + frac)
            perturbed = normalize_weights(perturbed)
            pert_cts = cts_with_weights(kinase_df, perturbed)
            rho, _ = spearmanr(baseline_cts.loc[kinases], pert_cts.loc[kinases])
            scenario_a_results.append({
                "perturbed_component": component, "fraction": frac,
                "spearman_rho": rho,
                "top10_overlap": top_k_overlap(baseline_cts, pert_cts, 10),
                "top20_overlap": top_k_overlap(baseline_cts, pert_cts, 20),
            })
    scenario_a_df = pd.DataFrame(scenario_a_results)
    scenario_a_df.to_csv(OUTPUT_DIR / "cts_weight_perturbation_scenarioA.csv", index=False)
    print("\n=== Scenario A: CTS one-at-a-time +/-25% perturbations ===")
    print(scenario_a_df.to_string(index=False))

    # ---------------- Scenario B: 1000 random weight vectors (CTS), constrained to sum=1 ----------------
    rng = np.random.default_rng(RNG_SEED)
    n_random = 1000
    components = list(CTS_WEIGHTS.keys())
    random_results = []
    for i in range(n_random):
        raw_w = rng.dirichlet(np.ones(len(components)) * 5)  # concentration=5: moderate spread around uniform
        w = dict(zip(components, raw_w))
        pert_cts = cts_with_weights(kinase_df, w)
        rho, _ = spearmanr(baseline_cts.loc[kinases], pert_cts.loc[kinases])
        random_results.append({
            "trial": i, "spearman_rho": rho,
            "top10_overlap": top_k_overlap(baseline_cts, pert_cts, 10),
            "top20_overlap": top_k_overlap(baseline_cts, pert_cts, 20),
            **{f"w_{c}": w[c] for c in components},
        })
    random_df = pd.DataFrame(random_results)
    random_df.to_csv(OUTPUT_DIR / "cts_weight_perturbation_random1000.csv", index=False)
    print(f"\n=== Scenario B: {n_random} random CTS weight draws ===")
    print(f"Spearman rho: mean={random_df['spearman_rho'].mean():.4f}  median={random_df['spearman_rho'].median():.4f}  min={random_df['spearman_rho'].min():.4f}")
    print(f"Top-10 overlap: mean={random_df['top10_overlap'].mean():.4f}  median={random_df['top10_overlap'].median():.4f}  min={random_df['top10_overlap'].min():.4f}")
    print(f"Top-20 overlap: mean={random_df['top20_overlap'].mean():.4f}  median={random_df['top20_overlap'].median():.4f}  min={random_df['top20_overlap'].min():.4f}")

    # ---------------- PairCTS weight perturbation: +/-25% one-at-a-time ----------------
    paircts_scenario_results = []
    for component in PAIR_CTS_WEIGHTS:
        for frac in (-0.25, 0.25):
            perturbed = dict(PAIR_CTS_WEIGHTS)
            perturbed[component] = PAIR_CTS_WEIGHTS[component] * (1 + frac)
            perturbed = normalize_weights(perturbed)
            pert_pairs = pairct_with_weights(baseline_cts, community_map, crosstalk_edges, perturbed, kinases)
            rho, _ = spearmanr(baseline_paircts, pert_pairs.loc[baseline_paircts.index])
            paircts_scenario_results.append({
                "perturbed_component": component, "fraction": frac,
                "spearman_rho": rho,
                "top10_overlap": top_k_overlap(baseline_paircts, pert_pairs, 10),
                "top20_overlap": top_k_overlap(baseline_paircts, pert_pairs, 20),
            })
    paircts_scenario_df = pd.DataFrame(paircts_scenario_results)
    paircts_scenario_df.to_csv(OUTPUT_DIR / "paircts_weight_perturbation_scenarioA.csv", index=False)
    print("\n=== PairCTS one-at-a-time +/-25% perturbations ===")
    print(paircts_scenario_df.to_string(index=False))

    # ---------------- PairCTS: 1000 random weight vectors ----------------
    pair_components = list(PAIR_CTS_WEIGHTS.keys())
    pair_random_results = []
    for i in range(n_random):
        raw_w = rng.dirichlet(np.ones(len(pair_components)) * 5)
        w = dict(zip(pair_components, raw_w))
        pert_pairs = pairct_with_weights(baseline_cts, community_map, crosstalk_edges, w, kinases)
        rho, _ = spearmanr(baseline_paircts, pert_pairs.loc[baseline_paircts.index])
        pair_random_results.append({
            "trial": i, "spearman_rho": rho,
            "top10_overlap": top_k_overlap(baseline_paircts, pert_pairs, 10),
            "top20_overlap": top_k_overlap(baseline_paircts, pert_pairs, 20),
        })
    pair_random_df = pd.DataFrame(pair_random_results)
    pair_random_df.to_csv(OUTPUT_DIR / "paircts_weight_perturbation_random1000.csv", index=False)
    print(f"\n=== PairCTS: {n_random} random weight draws ===")
    print(f"Spearman rho: mean={pair_random_df['spearman_rho'].mean():.4f}  median={pair_random_df['spearman_rho'].median():.4f}  min={pair_random_df['spearman_rho'].min():.4f}")
    print(f"Top-10 overlap: mean={pair_random_df['top10_overlap'].mean():.4f}  min={pair_random_df['top10_overlap'].min():.4f}")
    print(f"Top-20 overlap: mean={pair_random_df['top20_overlap'].mean():.4f}  min={pair_random_df['top20_overlap'].min():.4f}")

    # Save baseline rankings for downstream ablation/null-model scripts
    baseline_cts.to_csv(OUTPUT_DIR / "baseline_cts.csv", header=["CTS"])
    baseline_paircts.rename_axis(["kinase_i", "kinase_j"]).reset_index(name="PairCTS").to_csv(OUTPUT_DIR / "baseline_paircts.csv", index=False)


if __name__ == "__main__":
    main()
