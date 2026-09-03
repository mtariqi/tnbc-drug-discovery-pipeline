"""
Analysis 4: LMTK2/LMTK3 exclusion sensitivity check, and consistently-top-ranked
kinase/pair stability across the 1,000 random weight draws already run in
Analysis 1.

Two things this script does, both previously run only as one-off code and
never packaged into a reusable script:

1. LMTK2/LMTK3 exclusion: re-runs CTS and PairCTS on an 88-kinase panel
   (excluding LMTK2/LMTK3, per the kinase-panel audit in docs/limitations.md)
   and compares to the real 90-kinase baseline. Confirms whether excluding
   these two recently-reclassified genes changes any reported ranking.

2. Consistently-top-ranked stability: re-runs the SAME 1,000 random weight
   draws as Analysis 1 (identical RNG seed, so results are directly
   comparable), but this time tracks which specific kinases/pairs land in
   the top-20 on EACH trial, rather than only the aggregate correlation
   statistics Analysis 1 reports. Answers "which targets/pairs survive
   essentially every reasonable reweighting?" rather than just "how
   correlated is the ranking overall?"
"""
from __future__ import annotations

import itertools
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from kinase_scoring_pipeline import compute_cts, pair_cts, CTS_WEIGHTS, PAIR_CTS_WEIGHTS
from analysis1_weight_perturbation import load_real_inputs, top_k_overlap, RNG_SEED, BASE_DIR, OUTPUT_DIR

N_TRIALS = 1000
EXCLUDED_KINASES = {"LMTK2", "LMTK3"}


def run_lmtk_exclusion(kinase_df, crosstalk_edges, community_map, baseline_cts):
    kinases_full = list(kinase_df.index)
    kinases_reduced = [k for k in kinases_full if k not in EXCLUDED_KINASES]
    kinase_df_reduced = kinase_df.loc[kinases_reduced]

    reduced_cts = compute_cts(kinase_df_reduced, weights=CTS_WEIGHTS)["CTS"]
    rho, _ = spearmanr(baseline_cts.loc[kinases_reduced], reduced_cts.loc[kinases_reduced])
    t10 = top_k_overlap(baseline_cts.loc[kinases_reduced], reduced_cts, 10)
    t20 = top_k_overlap(baseline_cts.loc[kinases_reduced], reduced_cts, 20)

    baseline_top20 = set(baseline_cts.sort_values(ascending=False).head(20).index)
    lmtk_in_top20 = bool(EXCLUDED_KINASES & baseline_top20)

    def pairct_series(cts, kinases):
        rows = {}
        for ki, kj in itertools.combinations(kinases, 2):
            rows[(ki, kj)] = pair_cts(ki, kj, cts, community_map, crosstalk_edges)
        return pd.Series(rows)

    baseline_paircts_reduced_universe = pairct_series(baseline_cts, kinases_reduced)
    reduced_paircts = pairct_series(reduced_cts, kinases_reduced)
    rho_p, _ = spearmanr(baseline_paircts_reduced_universe, reduced_paircts.loc[baseline_paircts_reduced_universe.index])
    t10_p = top_k_overlap(baseline_paircts_reduced_universe, reduced_paircts, 10)
    t20_p = top_k_overlap(baseline_paircts_reduced_universe, reduced_paircts, 20)

    result = pd.DataFrame([
        {"analysis": "CTS_LMTK_exclusion", "spearman_rho": rho, "top10_overlap": t10, "top20_overlap": t20,
         "lmtk_in_baseline_top20": lmtk_in_top20},
        {"analysis": "PairCTS_LMTK_exclusion", "spearman_rho": rho_p, "top10_overlap": t10_p, "top20_overlap": t20_p,
         "lmtk_in_baseline_top20": None},
    ])
    return result


def run_stability_analysis(kinase_df, crosstalk_edges, community_map, baseline_cts):
    kinases = list(kinase_df.index)
    rng = np.random.default_rng(RNG_SEED)  # same seed as Analysis 1, for direct comparability

    baseline_top20_kinases = set(baseline_cts.sort_values(ascending=False).head(20).index)

    # ---- CTS: which kinases stay in top-20 across all N_TRIALS random draws ----
    components = list(CTS_WEIGHTS.keys())
    kinase_top20_counts = pd.Series(0, index=kinases)
    for i in range(N_TRIALS):
        w = dict(zip(components, rng.dirichlet(np.ones(len(components)) * 5)))
        cts_i = compute_cts(kinase_df, weights=w)["CTS"]
        top20 = set(cts_i.sort_values(ascending=False).head(20).index)
        for k in top20:
            kinase_top20_counts[k] += 1

    kinase_stability_df = kinase_top20_counts.sort_values(ascending=False).reset_index()
    kinase_stability_df.columns = ["kinase_id", "n_trials_in_top20_of_1000"]
    kinase_stability_df["in_baseline_top20"] = kinase_stability_df["kinase_id"].isin(baseline_top20_kinases)

    # ---- PairCTS: which pairs stay in top-20 across N_TRIALS random draws ----
    # Re-seed identically so the CTS-side and PairCTS-side draws match Analysis 1's
    # two separate random streams exactly.
    rng2 = np.random.default_rng(RNG_SEED)
    pair_components = list(PAIR_CTS_WEIGHTS.keys())
    all_pairs = list(itertools.combinations(kinases, 2))
    pair_top20_counts: dict = {}
    for i in range(N_TRIALS):
        w = dict(zip(pair_components, rng2.dirichlet(np.ones(len(pair_components)) * 5)))
        scores = {(ki, kj): pair_cts(ki, kj, baseline_cts, community_map, crosstalk_edges, weights=w) for ki, kj in all_pairs}
        s = pd.Series(scores)
        top20 = set(s.sort_values(ascending=False).head(20).index)
        for pair in top20:
            pair_top20_counts[pair] = pair_top20_counts.get(pair, 0) + 1

    pair_rows = [{"kinase_i": p[0], "kinase_j": p[1], "n_trials_in_top20_of_1000": c} for p, c in pair_top20_counts.items()]
    pair_stability_df = pd.DataFrame(pair_rows).sort_values("n_trials_in_top20_of_1000", ascending=False)

    return kinase_stability_df, pair_stability_df


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    kinase_df, crosstalk_edges, community_map = load_real_inputs()
    baseline_cts = compute_cts(kinase_df, weights=CTS_WEIGHTS)["CTS"]

    print("=== LMTK2/LMTK3 exclusion sensitivity check ===")
    lmtk_result = run_lmtk_exclusion(kinase_df, crosstalk_edges, community_map, baseline_cts)
    print(lmtk_result.to_string(index=False))
    lmtk_result.to_csv(OUTPUT_DIR / "lmtk_exclusion_analysis.csv", index=False)

    print("\n=== Consistently-top-ranked stability (1,000 random weight draws) ===")
    kinase_stability_df, pair_stability_df = run_stability_analysis(kinase_df, crosstalk_edges, community_map, baseline_cts)

    always_top20 = kinase_stability_df[kinase_stability_df["n_trials_in_top20_of_1000"] == N_TRIALS]
    print(f"\nCTS kinases in top-20 in ALL {N_TRIALS} trials ({len(always_top20)}):")
    print(always_top20["kinase_id"].tolist())
    print("\nTop 15 by stability:")
    print(kinase_stability_df.head(15).to_string(index=False))
    kinase_stability_df.to_csv(OUTPUT_DIR / "cts_kinase_top20_stability.csv", index=False)

    always_top20_pairs = pair_stability_df[pair_stability_df["n_trials_in_top20_of_1000"] == N_TRIALS]
    print(f"\nPairCTS pairs in top-20 in ALL {N_TRIALS} trials ({len(always_top20_pairs)}):")
    print(always_top20_pairs.to_string(index=False))
    print("\nTop 15 by stability:")
    print(pair_stability_df.head(15).to_string(index=False))
    pair_stability_df.to_csv(OUTPUT_DIR / "paircts_pair_top20_stability.csv", index=False)

    print(f"\nWrote lmtk_exclusion_analysis.csv, cts_kinase_top20_stability.csv, "
          f"paircts_pair_top20_stability.csv to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
