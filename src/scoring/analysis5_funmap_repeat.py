"""
Analysis 5: FunMap-substituted repeat of the PairCTS weight-perturbation,
ablation, and null-model analyses -- confirms the STRING-based findings in
Analyses 1-3 aren't an artifact of STRING's specific topology.

Design (per funmap_dual_network_scoring.py): CTS values and Louvain
community assignments are held fixed (from STRING, as elsewhere in this
project); only the crosstalk edge source is swapped from STRING to FunMap.

Requires funmap_90kinase_scored.tsv (the real, panel-restricted FunMap
network -- gene1/gene2/score columns) at
<repo>/data/raw/funmap/funmap_90kinase_scored.tsv. If your copy lives
elsewhere, set the FUNMAP_PATH environment variable.

Note on Null Model B: FunMap's real edge density within this panel is
~77%, too dense for degree-preserving double-edge-swap rewiring to be
tractable (few valid non-edges remain to swap into). This script instead
shuffles which real edge weight is assigned to which edge pair, holding
the near-complete topology fixed -- a necessary adaptation given the
network's density, documented in the manuscript (Section 4.11).
"""
from __future__ import annotations

import itertools
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from kinase_scoring_pipeline import compute_cts, pair_cts, CTS_WEIGHTS, PAIR_CTS_WEIGHTS
from funmap_dual_network_scoring import load_funmap_network
from analysis1_weight_perturbation import load_real_inputs, top_k_overlap, normalize_weights, RNG_SEED, BASE_DIR, OUTPUT_DIR
from analysis3_null_models import empirical_p_and_z

N_TRIALS = 1000
FUNMAP_PATH = Path(os.environ.get("FUNMAP_PATH", BASE_DIR / "data/raw/funmap/funmap_90kinase_scored.tsv"))


def pairct_series(kinases, cts, community_map, edges, weights=PAIR_CTS_WEIGHTS):
    rows = {}
    for ki, kj in itertools.combinations(kinases, 2):
        rows[(ki, kj)] = pair_cts(ki, kj, cts, community_map, edges, weights=weights)
    return pd.Series(rows)


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    kinase_df, string_edges, community_map = load_real_inputs()
    kinases = list(kinase_df.index)
    baseline_cts = compute_cts(kinase_df, weights=CTS_WEIGHTS)["CTS"]
    rng = np.random.default_rng(RNG_SEED)

    print(f"Loading real FunMap network from {FUNMAP_PATH}...")
    funmap_edges = load_funmap_network(str(FUNMAP_PATH), kinase_panel=kinases)
    print(f"Loaded {len(funmap_edges)} FunMap edges within the {len(kinases)}-kinase panel")

    string_baseline = pairct_series(kinases, baseline_cts, community_map, string_edges)
    funmap_baseline = pairct_series(kinases, baseline_cts, community_map, funmap_edges)

    rho, _ = spearmanr(string_baseline, funmap_baseline)
    t10 = top_k_overlap(string_baseline, funmap_baseline, 10)
    t20 = top_k_overlap(string_baseline, funmap_baseline, 20)
    print(f"\n=== STRING-based vs FunMap-based PairCTS rank stability ===")
    print(f"Spearman rho={rho:.4f}  top10_overlap={t10:.3f}  top20_overlap={t20:.3f}")

    # ---- Weight perturbation: 1,000 random draws, FunMap-based PairCTS ----
    pair_components = list(PAIR_CTS_WEIGHTS.keys())
    funmap_random_results = []
    for i in range(N_TRIALS):
        w = dict(zip(pair_components, rng.dirichlet(np.ones(len(pair_components)) * 5)))
        pert = pairct_series(kinases, baseline_cts, community_map, funmap_edges, weights=w)
        r, _ = spearmanr(funmap_baseline, pert.loc[funmap_baseline.index])
        funmap_random_results.append({
            "trial": i, "spearman_rho": r,
            "top10_overlap": top_k_overlap(funmap_baseline, pert, 10),
            "top20_overlap": top_k_overlap(funmap_baseline, pert, 20),
        })
        if (i + 1) % 250 == 0:
            print(f"  weight perturbation: {i+1}/{N_TRIALS}")
    funmap_random_df = pd.DataFrame(funmap_random_results)
    funmap_random_df.to_csv(OUTPUT_DIR / "paircts_funmap_weight_perturbation_random1000.csv", index=False)
    print(f"\nFunMap PairCTS weight perturbation: mean rho={funmap_random_df['spearman_rho'].mean():.4f} "
          f"median={funmap_random_df['spearman_rho'].median():.4f} min={funmap_random_df['spearman_rho'].min():.4f}")

    # ---- Ablation: remove complementarity / crosstalk, FunMap-based ----
    ablation_results = []
    for component in ("gamma", "delta"):
        ablated = {k: v for k, v in PAIR_CTS_WEIGHTS.items() if k != component}
        ablated = normalize_weights(ablated)
        full_ablated = {k: (0.0 if k == component else ablated[k]) for k in PAIR_CTS_WEIGHTS}
        ablated_scores = pairct_series(kinases, baseline_cts, community_map, funmap_edges, weights=full_ablated)
        r, _ = spearmanr(funmap_baseline, ablated_scores.loc[funmap_baseline.index])
        label = "complementarity" if component == "gamma" else "crosstalk_funmap"
        ablation_results.append({
            "ablated_component": label, "spearman_rho": r,
            "top10_overlap": top_k_overlap(funmap_baseline, ablated_scores, 10),
            "top20_overlap": top_k_overlap(funmap_baseline, ablated_scores, 20),
        })
    ablation_df = pd.DataFrame(ablation_results)
    ablation_df.to_csv(OUTPUT_DIR / "paircts_funmap_ablation_results.csv", index=False)
    print(f"\n=== FunMap-based PairCTS ablation ===")
    print(ablation_df.to_string(index=False))

    # ---- Null model: edge-weight shuffle (topology too dense to rewire) ----
    def best_pair_score(edges):
        best = -np.inf
        for ki, kj in itertools.combinations(kinases, 2):
            s = pair_cts(ki, kj, baseline_cts, community_map, edges)
            best = max(best, s)
        return best

    observed_funmap_best = best_pair_score(funmap_edges)
    print(f"\nObserved best FunMap-based PairCTS: {observed_funmap_best:.6f}")

    edge_keys = list(funmap_edges.keys())
    edge_weights = list(funmap_edges.values())
    null_scores = []
    for i in range(N_TRIALS):
        shuffled = list(edge_weights)
        rng.shuffle(shuffled)
        rand_edges = dict(zip(edge_keys, shuffled))
        null_scores.append(best_pair_score(rand_edges))
        if (i + 1) % 250 == 0:
            print(f"  null model: {i+1}/{N_TRIALS}")

    null_arr = np.array(null_scores)
    p_val, z_val = empirical_p_and_z(observed_funmap_best, null_arr)
    print(f"\n=== Null Model B' (FunMap: edge-weight shuffle, fixed near-complete topology) ===")
    print(f"Null mean={null_arr.mean():.4f} sd={null_arr.std():.4f}  Observed={observed_funmap_best:.4f}  "
          f"Z={z_val:.2f}  p={p_val:.4f}")
    pd.Series(null_arr, name="null_score").to_csv(OUTPUT_DIR / "null_model_B_funmap_weight_shuffle.csv", index=False)

    print(f"\nAll FunMap-repeat outputs written to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
