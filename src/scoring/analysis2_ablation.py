"""
Analysis 2: Leave-one-component-out ablation for CTS and PairCTS.

For each component, zero out its weight (redistributing nothing -- just set
to 0 and renormalize the rest) and compare the resulting ranking to baseline.
"""
from __future__ import annotations

import itertools
import pandas as pd
from scipy.stats import spearmanr

from kinase_scoring_pipeline import compute_cts, pair_cts, CTS_WEIGHTS, PAIR_CTS_WEIGHTS
from analysis1_weight_perturbation import load_real_inputs, top_k_overlap, normalize_weights, OUTPUT_DIR


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    kinase_df, crosstalk_edges, community_map = load_real_inputs()
    kinases = list(kinase_df.index)

    baseline_cts = compute_cts(kinase_df, weights=CTS_WEIGHTS)["CTS"]

    def pairct_series(cts, weights):
        rows = {}
        for ki, kj in itertools.combinations(kinases, 2):
            rows[(ki, kj)] = pair_cts(ki, kj, cts, community_map, crosstalk_edges, weights=weights)
        return pd.Series(rows)

    baseline_paircts = pairct_series(baseline_cts, PAIR_CTS_WEIGHTS)

    # ---- CTS ablations: drop each of the four components ----
    cts_results = []
    for component in CTS_WEIGHTS:
        ablated = {k: v for k, v in CTS_WEIGHTS.items() if k != component}
        ablated = normalize_weights(ablated)
        # compute_cts always computes all four terms internally; to "remove" one,
        # set its weight to 0 explicitly (not omit the key, since compute_cts
        # indexes CTS_WEIGHTS by name and expects all four keys present).
        full_ablated_weights = {k: (0.0 if k == component else ablated[k]) for k in CTS_WEIGHTS}
        ablated_cts = compute_cts(kinase_df, weights=full_ablated_weights)["CTS"]
        rho, _ = spearmanr(baseline_cts.loc[kinases], ablated_cts.loc[kinases])
        cts_results.append({
            "ablated_component": component,
            "spearman_rho": rho,
            "top10_overlap": top_k_overlap(baseline_cts, ablated_cts, 10),
            "top20_overlap": top_k_overlap(baseline_cts, ablated_cts, 20),
        })
    cts_ablation_df = pd.DataFrame(cts_results).sort_values("spearman_rho")
    cts_ablation_df.to_csv(OUTPUT_DIR / "cts_ablation_results.csv", index=False)
    print("=== CTS leave-one-out ablation (sorted by impact, most disruptive first) ===")
    print(cts_ablation_df.to_string(index=False))

    # ---- PairCTS ablations: drop complementarity (gamma) or crosstalk (delta) ----
    pair_results = []
    for component in ("gamma", "delta"):
        ablated = {k: v for k, v in PAIR_CTS_WEIGHTS.items() if k != component}
        ablated = normalize_weights(ablated)
        full_ablated_weights = {k: (0.0 if k == component else ablated[k]) for k in PAIR_CTS_WEIGHTS}
        ablated_pairs = pairct_series(baseline_cts, full_ablated_weights)
        rho, _ = spearmanr(baseline_paircts, ablated_pairs.loc[baseline_paircts.index])
        label = "complementarity" if component == "gamma" else "crosstalk"
        pair_results.append({
            "ablated_component": label,
            "spearman_rho": rho,
            "top10_overlap": top_k_overlap(baseline_paircts, ablated_pairs, 10),
            "top20_overlap": top_k_overlap(baseline_paircts, ablated_pairs, 20),
        })
    # Also ablate alpha/beta (the CTS(i)/CTS(j) terms) for completeness, since
    # Dr Miah's spec focuses on complementarity/crosstalk but a full six-component
    # picture (matching "no essentiality...no crosstalk" framing) is more informative.
    for component, label in (("alpha", "target_cts_i"), ("beta", "target_cts_j")):
        ablated = {k: v for k, v in PAIR_CTS_WEIGHTS.items() if k != component}
        ablated = normalize_weights(ablated)
        full_ablated_weights = {k: (0.0 if k == component else ablated[k]) for k in PAIR_CTS_WEIGHTS}
        ablated_pairs = pairct_series(baseline_cts, full_ablated_weights)
        rho, _ = spearmanr(baseline_paircts, ablated_pairs.loc[baseline_paircts.index])
        pair_results.append({
            "ablated_component": label,
            "spearman_rho": rho,
            "top10_overlap": top_k_overlap(baseline_paircts, ablated_pairs, 10),
            "top20_overlap": top_k_overlap(baseline_paircts, ablated_pairs, 20),
        })
    pair_ablation_df = pd.DataFrame(pair_results).sort_values("spearman_rho")
    pair_ablation_df.to_csv(OUTPUT_DIR / "paircts_ablation_results.csv", index=False)
    print("\n=== PairCTS leave-one-out ablation (sorted by impact, most disruptive first) ===")
    print(pair_ablation_df.to_string(index=False))


if __name__ == "__main__":
    main()
