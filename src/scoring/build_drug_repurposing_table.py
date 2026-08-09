"""
build_drug_repurposing_table.py

Builds the "Drug Repurposing Opportunities" table for the manuscript's
top-ranked CTS kinases, using ONLY real DGIdb interaction data already
loaded elsewhere in this project (CTS_all_90_kinases.tsv, dgidb_interactions.tsv).
No scores, drugs, or approval statuses are invented -- every row in the
output table comes directly from a real row in the real DGIdb file.

This does not introduce a second scoring system: it is a direct lookup
against the CTS ranking (Section 4.1) and DGIdb (Section 2.1), not a new
prioritization method.

Usage (run from src/scoring/, alongside kinase_scoring_pipeline.py):
    python3 build_drug_repurposing_table.py [base_dir] [--top-n N]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


def load_top_cts_kinases(base_dir: str, top_n: int) -> pd.Series:
    """Real CTS ranking, already computed and stored (Section 4.1)."""
    cts_df = pd.read_csv(Path(base_dir) / "data/processed/CTS_all_90_kinases.tsv", sep="\t", index_col="kinase_id")
    cts = cts_df["CTS"].sort_values(ascending=False)
    return cts.head(top_n)


def load_real_dgidb_interactions(base_dir: str) -> pd.DataFrame:
    """Real DGIdb interaction records already downloaded for this project (Section 2.1)."""
    return pd.read_csv(Path(base_dir) / "data/processed/dgidb/dgidb_interactions.tsv", sep="\t")


def build_repurposing_table(
    top_kinases: pd.Series,
    dgidb_df: pd.DataFrame,
    max_drugs_per_kinase: int = 3,
) -> pd.DataFrame:
    """
    For each top-ranked kinase, lists up to max_drugs_per_kinase real drugs
    from the real DGIdb records for that kinase, ordered by source count
    (more independently-curated sources first) -- the same tie-breaking
    rule already used elsewhere in this project (build_real_gene_drugs.py),
    not a new ranking invented for this table.

    A kinase with zero real DGIdb-recorded drugs is listed with an
    explicit "No DGIdb-recorded inhibitor" note rather than omitted or
    filled in -- consistent with this project's standing rule of never
    silently hiding a real gap.
    """
    rows = []
    for kinase, cts_score in top_kinases.items():
        matches = dgidb_df[dgidb_df["kinase_id"] == kinase].copy()
        if matches.empty:
            rows.append({
                "Kinase": kinase, "CTS": round(cts_score, 3),
                "Candidate drug(s)": "No DGIdb-recorded inhibitor", "Source": "\u2014",
            })
            continue

        matches["n_sources"] = matches["sources"].fillna("").apply(
            lambda s: len([x for x in str(s).split(",") if x.strip()])
        )
        top_drugs = (
            matches.sort_values("n_sources", ascending=False)
            .drop_duplicates(subset="drug")
            .head(max_drugs_per_kinase)
        )
        rows.append({
            "Kinase": kinase,
            "CTS": round(cts_score, 3),
            "Candidate drug(s)": ", ".join(top_drugs["drug"].tolist()),
            "Source": "DGIdb",
        })
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("base_dir", nargs="?", default=str(Path.home() / "rtk_nrtk_tnbc"))
    ap.add_argument("--top-n", type=int, default=5)
    args = ap.parse_args()

    print(f"Loading real CTS ranking (top {args.top_n})...")
    top_kinases = load_top_cts_kinases(args.base_dir, args.top_n)
    print(top_kinases.to_string())

    print("\nLoading real DGIdb interactions...")
    dgidb_df = load_real_dgidb_interactions(args.base_dir)
    print(f"  {len(dgidb_df)} real interaction records loaded")

    print("\nBuilding drug repurposing table from real DGIdb records only...")
    table = build_repurposing_table(top_kinases, dgidb_df)
    print("\n" + table.to_string(index=False))

    out_path = Path(args.base_dir) / "data/processed/drug_repurposing_table.tsv"
    table.to_csv(out_path, sep="\t", index=False)
    print(f"\nWrote {out_path}")


# =====================================================================
# SMOKE TEST -- synthetic data, verifies real-data-only behavior:
# a kinase with real drugs gets them in the right order; a kinase with
# NO real DGIdb record is flagged explicitly, never silently invented.
# =====================================================================

def _run_smoke_test():
    top_kinases = pd.Series({"ERBB2": 0.690, "EGFR": 0.613, "PTK2": 0.532, "OBSCURE_KINASE": 0.410})

    dgidb_df = pd.DataFrame([
        {"kinase_id": "ERBB2", "drug": "trastuzumab", "sources": "DrugBank,ChEMBL,FDA"},
        {"kinase_id": "ERBB2", "drug": "pertuzumab", "sources": "DrugBank,FDA"},
        {"kinase_id": "ERBB2", "drug": "lapatinib", "sources": "DrugBank"},
        {"kinase_id": "ERBB2", "drug": "neratinib", "sources": "ChEMBL"},
        {"kinase_id": "EGFR", "drug": "erlotinib", "sources": "DrugBank,ChEMBL"},
        {"kinase_id": "EGFR", "drug": "gefitinib", "sources": "DrugBank"},
        {"kinase_id": "PTK2", "drug": "defactinib", "sources": "ChEMBL"},
        # OBSCURE_KINASE has no real record at all
    ])

    print("=== Testing real-data-only lookup (no invented drugs, no invented scores) ===")
    table = build_repurposing_table(top_kinases, dgidb_df, max_drugs_per_kinase=3)
    print(table.to_string(index=False))

    erbb2_row = table[table["Kinase"] == "ERBB2"].iloc[0]
    assert erbb2_row["CTS"] == 0.690
    # Should cap at 3 drugs, prioritizing by real source count: trastuzumab(3) > pertuzumab(2) > lapatinib/neratinib(1, tie)
    assert erbb2_row["Candidate drug(s)"].split(", ")[0] == "trastuzumab"
    assert len(erbb2_row["Candidate drug(s)"].split(", ")) == 3
    print("PASSED: ERBB2 correctly capped at 3 drugs, ordered by real source count, trastuzumab first.\n")

    obscure_row = table[table["Kinase"] == "OBSCURE_KINASE"].iloc[0]
    assert obscure_row["Candidate drug(s)"] == "No DGIdb-recorded inhibitor"
    print("PASSED: kinase with zero real DGIdb records is explicitly flagged, not silently omitted or fabricated.\n")

    print("=" * 70)
    print("ALL SMOKE TESTS PASSED. Ready to run against real CTS + DGIdb data.")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--test":
        _run_smoke_test()
    else:
        main()
