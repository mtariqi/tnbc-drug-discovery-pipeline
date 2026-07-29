"""
build_and_verify_real_cts.py

Assembles the real kinase_df compute_cts() expects, directly from your raw
data files (not the already-computed CTS_all_90_kinases.tsv), then compares
the freshly-computed result against that existing file as a real ground-truth
check -- if they match closely, this confirms the assembly logic reconstructs
the original methodology; if not, that's real signal something differs, not
something to paper over.

Honest gaps, stated up front:
    - chembl_count and trial_stage are NOT available in this data (no live
      ChEMBL API access from here, and no trial-stage source was ever built).
      Only dgidb_score is used for the druggability term -- compute_cts()'s
      "average of available drug_cols" degrades gracefully to just this one
      column rather than fabricating the other two.
    - dgidb_score is computed as the MEAN interaction_score per kinase, across
      only DGIdb rows that pass the already-validated junk-name filter (so
      "antineoplastic agent"-style noise doesn't inflate a kinase's apparent
      druggability). Mean (not sum) was chosen because interaction_score looks
      like a per-interaction confidence/strength value, not a count -- but
      this is a real methodological choice, not a certainty, which is exactly
      why this script compares against the existing real CTS_all_90_kinases.tsv
      rather than trusting the choice blindly.
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
# string_network_builder.py lives in the sibling src/data_loaders/ folder in this repo's
# structure, not next to this script in src/scoring/ -- add that path too.
sys.path.insert(0, str(Path(__file__).parent.parent / "data_loaders"))

from kinase_scoring_pipeline import compute_cts
from string_network_builder import compute_centrality

# Reuse the already-validated junk-name filter from run_pairs_and_triplets.py
try:
    from run_pairs_and_triplets import _is_junk_drug_name, _normalize_drug_name
except ImportError:
    raise SystemExit(
        "run_pairs_and_triplets.py (with the drug-name patch already applied) must be "
        "in the same directory -- copy it here first."
    )


def build_kinase_df(base_dir: str) -> pd.DataFrame:
    base = Path(base_dir)

    kinase_list = [k.strip() for k in open(base / "data/raw/kinases/kinase_90_list.txt") if k.strip()]
    df = pd.DataFrame(index=pd.Index(kinase_list, name="kinase_id"))

    # --- Essentiality ---
    depmap = pd.read_csv(base / "data/processed/depmap/depmap_tnbc_essentiality.tsv", sep="\t")
    depmap = depmap.set_index("kinase_id")
    df["depmap_score"] = depmap["mean_depmap_score"]

    # --- Survival ---
    survival = pd.read_csv(base / "data/processed/tcga_brca/survival_stats.tsv", sep="\t")
    survival = survival.set_index("kinase_id")
    df["cox_hr"] = survival["cox_hr"]
    df["logrank_p"] = survival["logrank_p"]

    # --- Centrality (real, already-validated function, not reimplemented) ---
    edges_df = pd.read_csv(base / "data/processed/string/string_edges.tsv", sep="\t")
    centrality_df = compute_centrality(edges_df, all_kinases=kinase_list)
    centrality_df = centrality_df.set_index("kinase_id")
    df["betweenness"] = centrality_df["betweenness"]
    df["pagerank"] = centrality_df["pagerank"]
    df["degree"] = centrality_df["degree"]

    # --- Druggability: mean interaction_score of non-junk drugs per kinase ---
    dgidb = pd.read_csv(base / "data/processed/dgidb/dgidb_interactions.tsv", sep="\t")
    dgidb_clean = dgidb[~dgidb["drug"].apply(_is_junk_drug_name)].copy()
    n_dropped = len(dgidb) - len(dgidb_clean)
    print(f"  Dropped {n_dropped}/{len(dgidb)} DGIdb rows as junk drug names before aggregating "
          f"druggability scores.")
    dgidb_score_mean = dgidb_clean.groupby("kinase_id")["interaction_score"].mean()
    dgidb_score_count = dgidb_clean.groupby("kinase_id")["drug"].nunique()
    df["dgidb_score"] = dgidb_score_mean          # original hypothesis: mean confidence
    df["dgidb_score_count"] = dgidb_score_count   # alternative hypothesis: breadth of evidence
    # chembl_count, trial_stage: genuinely not available -- not added, so compute_cts()
    # falls back to using dgidb_score alone for the druggability term.

    return df


def main():
    base_dir = sys.argv[1] if len(sys.argv) > 1 else str(Path.home() / "rtk_nrtk_tnbc")

    print("Building kinase_df from real raw data files...")
    kinase_df = build_kinase_df(base_dir)
    print(f"  Built kinase_df: {len(kinase_df)} kinases")

    existing_path = Path(base_dir) / "data/processed/CTS_all_90_kinases.tsv"
    existing = None
    if existing_path.exists():
        existing = pd.read_csv(existing_path, sep="\t", index_col="kinase_id")

    def run_variant(label, druggability_col):
        variant_df = kinase_df.drop(columns=["dgidb_score", "dgidb_score_count"]).copy()
        variant_df["dgidb_score"] = kinase_df[druggability_col]
        result = compute_cts(variant_df)
        print(f"\n=== Variant: {label} ===")
        print("Top 10 by fresh CTS:")
        print(result.sort_values("CTS", ascending=False)[["CTS"]].head(10).to_string())
        if existing is not None:
            compare = pd.DataFrame({"CTS_fresh": result["CTS"], "CTS_existing": existing["CTS"]})
            compare["abs_diff"] = (compare["CTS_fresh"] - compare["CTS_existing"]).abs()
            print(f"Mean abs diff: {compare['abs_diff'].mean():.4f}  |  Max abs diff: {compare['abs_diff'].max():.4f}")
        return result

    result_mean = run_variant("mean(interaction_score) -- original hypothesis", "dgidb_score")
    result_count = run_variant("count(distinct drugs) -- alternative hypothesis", "dgidb_score_count")

    if existing is not None:
        print("\n=== Which variant is closer, per-kinase? ===")
        compare_both = pd.DataFrame({
            "CTS_existing": existing["CTS"],
            "CTS_mean_variant": result_mean["CTS"],
            "CTS_count_variant": result_count["CTS"],
        })
        compare_both["mean_diff"] = (compare_both["CTS_mean_variant"] - compare_both["CTS_existing"]).abs()
        compare_both["count_diff"] = (compare_both["CTS_count_variant"] - compare_both["CTS_existing"]).abs()
        compare_both["count_variant_is_closer"] = compare_both["count_diff"] < compare_both["mean_diff"]
        print(f"count-variant closer for {compare_both['count_variant_is_closer'].sum()}/{len(compare_both)} kinases")
        print(f"mean-variant total abs diff:  {compare_both['mean_diff'].sum():.3f}")
        print(f"count-variant total abs diff: {compare_both['count_diff'].sum():.3f}")

    out_path = Path(base_dir) / "data/processed/CTS_all_90_kinases_FRESH.tsv"
    result_count.to_csv(out_path, sep="\t")
    print(f"\nWrote count-variant result to {out_path} (does not overwrite the existing file).")


if __name__ == "__main__":
    main()
