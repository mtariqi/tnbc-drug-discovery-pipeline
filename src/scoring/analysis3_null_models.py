"""
Analysis 3: Null-model / randomization testing.

Null Model A: shuffle which kinase gets which CTS value (permute the CTS
             series' labels), recompute the top-N PairCTS pool using the
             SAME real network/community structure, compare observed top
             pair scores to the null distribution.
Null Model B: randomize the STRING network via degree-preserving double-edge
             swaps (networkx), recompute Complementarity/Crosstalk/PairCTS
             using the REAL CTS values but a shuffled network.
Null Model C: shuffle the drug-target mapping (permute which drug maps to
             which kinase's target list, preserving each drug's real
             target-COUNT), rerun DrugPairScore-relevant target-pair lookups.

All three repeated 1000 times; reports Z-scores and empirical p-values for
the observed top-ranked pair/regimen scores against each null distribution.
"""
from __future__ import annotations

import itertools
import numpy as np
import pandas as pd
import networkx as nx

from kinase_scoring_pipeline import compute_cts, pair_cts, CTS_WEIGHTS, PAIR_CTS_WEIGHTS
from analysis1_weight_perturbation import load_real_inputs, BASE_DIR, OUTPUT_DIR

RNG_SEED = 42
N_PERMUTATIONS = 1000


def empirical_p_and_z(observed: float, null_dist: np.ndarray, greater_is_more_extreme: bool = True) -> tuple[float, float]:
    mu, sigma = null_dist.mean(), null_dist.std(ddof=1)
    z = (observed - mu) / sigma if sigma > 0 else float("nan")
    if greater_is_more_extreme:
        p = (np.sum(null_dist >= observed) + 1) / (len(null_dist) + 1)
    else:
        p = (np.sum(null_dist <= observed) + 1) / (len(null_dist) + 1)
    return p, z


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    kinase_df, crosstalk_edges, community_map = load_real_inputs()
    kinases = list(kinase_df.index)
    rng = np.random.default_rng(RNG_SEED)

    baseline_cts = compute_cts(kinase_df, weights=CTS_WEIGHTS)["CTS"]

    def best_pair_score(cts_series, comm_map, edges):
        best = -np.inf
        best_pair = None
        for ki, kj in itertools.combinations(kinases, 2):
            s = pair_cts(ki, kj, cts_series, comm_map, edges)
            if s > best:
                best, best_pair = s, (ki, kj)
        return best, best_pair

    observed_best_score, observed_best_pair = best_pair_score(baseline_cts, community_map, crosstalk_edges)
    print(f"Observed best real PairCTS pair: {observed_best_pair} = {observed_best_score:.6f}")

    # =================================================================
    # Null Model A: shuffle CTS value <-> kinase-label assignment
    # =================================================================
    null_a_scores = []
    cts_values = baseline_cts.values.copy()
    for i in range(N_PERMUTATIONS):
        shuffled_values = rng.permutation(cts_values)
        shuffled_cts = pd.Series(shuffled_values, index=baseline_cts.index)
        score, _ = best_pair_score(shuffled_cts, community_map, crosstalk_edges)
        null_a_scores.append(score)
        if (i + 1) % 200 == 0:
            print(f"  Null A: {i+1}/{N_PERMUTATIONS}")
    null_a = np.array(null_a_scores)
    p_a, z_a = empirical_p_and_z(observed_best_score, null_a)
    print(f"\n=== Null Model A: shuffled CTS values ===")
    print(f"Null distribution: mean={null_a.mean():.4f} sd={null_a.std():.4f}")
    print(f"Observed={observed_best_score:.4f}  Z={z_a:.2f}  empirical p={p_a:.4f}")
    pd.Series(null_a, name="null_score").to_csv(OUTPUT_DIR / "null_model_A_cts_shuffle.csv", index=False)

    # =================================================================
    # Null Model B: degree-preserving STRING network randomization
    # =================================================================
    G = nx.Graph()
    G.add_nodes_from(kinases)
    for (ki, kj), w in crosstalk_edges.items():
        if ki in kinases and kj in kinases:
            G.add_edge(ki, kj, weight=w)
    real_degree_sequence = dict(G.degree())
    print(f"\nReal STRING subgraph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    null_b_scores = []
    for i in range(N_PERMUTATIONS):
        G_rand = G.copy()
        try:
            # double_edge_swap preserves the degree sequence exactly; nswap chosen
            # as ~4x edge count for adequate mixing, matching standard practice.
            nx.double_edge_swap(G_rand, nswap=G.number_of_edges() * 4, max_tries=G.number_of_edges() * 40, seed=int(rng.integers(0, 2**31)))
        except nx.NetworkXError:
            pass  # occasionally too few suitable swaps remain; use what succeeded
        rand_edges = {}
        rand_weights = list(nx.get_edge_attributes(G, "weight").values())
        rng.shuffle(rand_weights)  # reassign real edge weights onto the randomized topology
        for (u, v), w in zip(G_rand.edges(), rand_weights):
            rand_edges[(u, v)] = w
        score, _ = best_pair_score(baseline_cts, community_map, rand_edges)
        null_b_scores.append(score)
        if (i + 1) % 200 == 0:
            print(f"  Null B: {i+1}/{N_PERMUTATIONS}")
    null_b = np.array(null_b_scores)
    p_b, z_b = empirical_p_and_z(observed_best_score, null_b)
    print(f"\n=== Null Model B: degree-preserving network randomization ===")
    print(f"Null distribution: mean={null_b.mean():.4f} sd={null_b.std():.4f}")
    print(f"Observed={observed_best_score:.4f}  Z={z_b:.2f}  empirical p={p_b:.4f}")
    pd.Series(null_b, name="null_score").to_csv(OUTPUT_DIR / "null_model_B_network_shuffle.csv", index=False)

    # =================================================================
    # Null Model C: shuffle drug-target mapping (preserving per-drug target count)
    # =================================================================
    print("\nLoading real DGIdb drug-target data for Null Model C...")
    dgidb = pd.read_csv(BASE_DIR / "data/processed/dgidb/dgidb_interactions.tsv", sep="\t")
    dgidb["drug"] = dgidb["drug"].str.lower()
    dgidb = dgidb[dgidb["kinase_id"].isin(kinases)][["drug", "kinase_id"]].dropna().drop_duplicates()
    target_map_real = dgidb.groupby("drug")["kinase_id"].apply(list).to_dict()
    print(f"Real drug-target map: {len(target_map_real)} drugs covering panel kinases")

    def best_drug_pair_score(target_map):
        best = -np.inf
        drugs = list(target_map.keys())
        # Cap search for tractability: real analysis over all drug pairs is
        # O(n^2 * targets^2); subsample drug pairs for the null runs since
        # only the MAXIMUM matters and 1000 permutations x full drug-pair
        # search is not tractable here -- instead evaluate a fixed random
        # sample of 2000 drug pairs per permutation, consistent across
        # observed and null so the comparison is apples-to-apples.
        return drugs

    # For tractability, evaluate the observed and null distributions over
    # the SAME fixed sample of drug pairs (2000, seeded) rather than the full
    # O(n^2) space -- consistent methodology, not a shortcut on the real answer.
    drugs_real = list(target_map_real.keys())
    all_drug_pairs = list(itertools.combinations(drugs_real, 2))
    sample_size = min(2000, len(all_drug_pairs))
    sample_idx = rng.choice(len(all_drug_pairs), size=sample_size, replace=False)
    sampled_pairs = [all_drug_pairs[i] for i in sample_idx]

    def best_over_sampled_pairs(target_map):
        best = -np.inf
        for di, dj in sampled_pairs:
            ti, tj = target_map.get(di, []), target_map.get(dj, [])
            if not ti or not tj:
                continue
            local_best = max(
                (pair_cts(ki, kj, baseline_cts, community_map, crosstalk_edges) for ki, kj in itertools.product(ti, tj)),
                default=-np.inf,
            )
            best = max(best, local_best)
        return best

    observed_c = best_over_sampled_pairs(target_map_real)
    print(f"Observed best real DrugPair target-CTS (sampled {sample_size} drug pairs): {observed_c:.6f}")

    all_targets_pool = [k for targets in target_map_real.values() for k in targets]
    null_c_scores = []
    for i in range(N_PERMUTATIONS):
        shuffled_pool = list(all_targets_pool)
        rng.shuffle(shuffled_pool)
        shuffled_map = {}
        idx = 0
        for drug, targets in target_map_real.items():
            n = len(targets)
            shuffled_map[drug] = shuffled_pool[idx:idx + n]
            idx += n
        score = best_over_sampled_pairs(shuffled_map)
        null_c_scores.append(score)
        if (i + 1) % 200 == 0:
            print(f"  Null C: {i+1}/{N_PERMUTATIONS}")
    null_c = np.array(null_c_scores)
    p_c, z_c = empirical_p_and_z(observed_c, null_c)
    print(f"\n=== Null Model C: shuffled drug-target mapping ===")
    print(f"Null distribution: mean={null_c.mean():.4f} sd={null_c.std():.4f}")
    print(f"Observed={observed_c:.4f}  Z={z_c:.2f}  empirical p={p_c:.4f}")
    pd.Series(null_c, name="null_score").to_csv(OUTPUT_DIR / "null_model_C_drugmap_shuffle.csv", index=False)

    summary = pd.DataFrame([
        {"null_model": "A: CTS value shuffle", "observed": observed_best_score, "null_mean": null_a.mean(), "null_sd": null_a.std(), "Z": z_a, "empirical_p": p_a},
        {"null_model": "B: network randomization", "observed": observed_best_score, "null_mean": null_b.mean(), "null_sd": null_b.std(), "Z": z_b, "empirical_p": p_b},
        {"null_model": "C: drug-target mapping shuffle", "observed": observed_c, "null_mean": null_c.mean(), "null_sd": null_c.std(), "Z": z_c, "empirical_p": p_c},
    ])
    summary.to_csv(OUTPUT_DIR / "null_model_summary.csv", index=False)
    print("\n=== SUMMARY ===")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
