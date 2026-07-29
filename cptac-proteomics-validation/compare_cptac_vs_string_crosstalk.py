"""
compare_cptac_vs_string_crosstalk.py

Compares REAL, MEASURED protein co-expression correlations (CPTAC BRCA
proteomics, this repo's cptac-proteomics-validation/) against the PREDICTED
kinase-kinase crosstalk edges from STRING already used in PairCTS/TripletCTS
scoring (src/scoring/kinase_scoring_pipeline.py's crosstalk_strength()).

WHY THIS MATTERS: everywhere else in this project, "crosstalk" between two
kinases is a STRING-predicted network property (co-occurrence in curated
databases, text-mining, etc.) -- never checked against real, measured
co-expression data. This is the one place in the whole project that can
make that comparison directly.

HONEST LIMITATION, stated up front (not new -- already documented in
docs/limitations.md): this CPTAC cohort is NOT restricted to TNBC -- it's
BRCA-wide (151 samples, all subtypes mixed). A real correlation here
reflects general breast cancer co-regulation, not necessarily TNBC-specific
biology. Treat agreement/disagreement with STRING as suggestive, not as a
TNBC-specific validation.

Run from cptac-proteomics-validation/ (needs src/data_loaders/'s real
string_edges.tsv path):
    python3 compare_cptac_vs_string_crosstalk.py
"""

import sys
from pathlib import Path

import pandas as pd


def extract_gene_symbol(compound_id: str) -> str:
    """
    CPTAC's raw gene identifier is a pipe-delimited compound string
    (ENSP...|ENST...|ENSG...|OTTHUMG...|OTTHUMT...|SYMBOL-201|SYMBOL|length).
    The real gene symbol is the second-to-last field -- confirmed against
    real data (e.g. 'BTK', 'ERBB3', 'ABL1' all extracted correctly), not
    the raw string correlation_analysis_final.py's text summary mistakenly
    printed (a cosmetic bug in that one script's print statements only --
    its own CSV outputs, and this extraction, are correct).
    """
    parts = compound_id.split("|")
    return parts[-2] if len(parts) >= 2 else compound_id


def main():
    base_dir = sys.argv[1] if len(sys.argv) > 1 else str(Path.home() / "rtk_nrtk_tnbc")
    string_edges_path = Path(base_dir) / "data/processed/string/string_edges.tsv"

    if not string_edges_path.exists():
        raise SystemExit(f"{string_edges_path} not found -- adjust base_dir or the path above.")

    print(f"Loading real STRING crosstalk edges from {string_edges_path}...")
    string_edges = pd.read_csv(string_edges_path, sep="\t")
    string_pairs = set()
    for _, row in string_edges.iterrows():
        string_pairs.add(tuple(sorted([row["source"], row["target"]])))
    print(f"  {len(string_pairs)} STRING-predicted kinase-kinase edges")

    cptac_path = Path(__file__).parent / "results" / "significant_correlations.csv"
    print(f"\nLoading real CPTAC significant correlations from {cptac_path}...")
    cptac = pd.read_csv(cptac_path)
    cptac["gene1"] = cptac["Gene1"].apply(extract_gene_symbol)
    cptac["gene2"] = cptac["Gene2"].apply(extract_gene_symbol)
    cptac["pair"] = cptac.apply(lambda r: tuple(sorted([r["gene1"], r["gene2"]])), axis=1)
    print(f"  {len(cptac)} statistically significant CPTAC correlation pairs")

    cptac["in_string"] = cptac["pair"].apply(lambda p: p in string_pairs)

    print(f"\n=== Overlap ===")
    n_in_string = cptac["in_string"].sum()
    print(f"CPTAC-significant pairs that ALSO have a STRING edge: {n_in_string}/{len(cptac)}")

    print(f"\n=== Real measured correlations that DO have a matching STRING edge "
          f"(real crosstalk_strength() would apply to these) ===")
    matched = cptac[cptac["in_string"]].sort_values("Correlation", key=abs, ascending=False)
    print(matched[["gene1", "gene2", "Correlation", "FDR"]].head(15).to_string(index=False))

    print(f"\n=== Real measured correlations with NO matching STRING edge "
          f"(crosstalk_strength() currently credits these as 0, despite real co-expression) ===")
    unmatched = cptac[~cptac["in_string"]].sort_values("Correlation", key=abs, ascending=False)
    print(unmatched[["gene1", "gene2", "Correlation", "FDR"]].head(15).to_string(index=False))

    out_path = Path(__file__).parent / "results" / "cptac_vs_string_comparison.csv"
    cptac[["gene1", "gene2", "Correlation", "FDR", "in_string"]].to_csv(out_path, index=False)
    print(f"\nWrote full comparison to {out_path}")

    print("\nHONEST REMINDER: this CPTAC cohort is BRCA-wide, not TNBC-restricted (see "
          "docs/limitations.md). A real correlation here is suggestive of general breast "
          "cancer co-regulation, not confirmed TNBC-specific biology.")


if __name__ == "__main__":
    main()
