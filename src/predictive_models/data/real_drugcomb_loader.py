"""
Real DrugComb loader -- UNVERIFIED SCHEMA, inspect before trusting.

DrugComb's public bulk download (https://drugcomb.org/download/) has historically
shipped a "summary" CSV with columns resembling: block_id, drug_row, drug_col,
cell_line_name, tissue_name, study_name, synergy_zip, synergy_bliss, synergy_loewe,
synergy_hsa, ic50_row, ic50_col, css_row, css_col, S (an overall score). This module
is written against that expected schema, but -- following the same discipline as
depmap_multiomic_loader.py's inspect_depmap_csv() -- it has NOT been confirmed
against your actual downloaded file this session. Run inspect_drugcomb_csv() first.

WHAT THIS DOES NOT DO: fetch or match SMILES strings for each drug. DrugComb's
summary file identifies drugs by name/ID, not by structure -- you'll need a
separate drug-name -> SMILES mapping (e.g. from DrugComb's own drug metadata file,
or PubChem/ChEMBL lookup) before src.predictive_models.encoders.molecule.featurize_smiles can
run. That join is also unverified and not implemented here -- do it as an explicit,
inspectable step rather than silently inside this loader.

WHAT THIS DOES NOT HANDLE: biologics (e.g. trastuzumab) with no small-molecule
SMILES. See models/synergy_gnn.py's docstring for why, and decide how to handle
biologic-containing rows (drop, or a fallback ID-embedding path) before training.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def inspect_drugcomb_csv(path: str, n_rows: int = 5) -> None:
    """Print columns and a sample so you can confirm this loader's column-name
    assumptions before trusting it -- same pattern as inspect_depmap_csv()."""
    df = pd.read_csv(path, nrows=n_rows)
    print(f"=== {path} ===")
    print(f"Columns: {df.columns.tolist()}")
    print(df.head(n_rows))


def load_drugcomb_summary(
    csv_path: str,
    drug_name_to_smiles: dict[str, str],
    cell_lines: list[str] | None = None,
    synergy_columns: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Load and restrict DrugComb's summary file to rows where BOTH drugs have a
    known SMILES mapping (required for featurization) and, optionally, a specific
    cell-line subset.

    `synergy_columns` maps your desired internal metric name -> the real column
    name in the CSV (default guesses below -- CONFIRM against inspect_drugcomb_csv()
    output first).
    """
    synergy_columns = synergy_columns or {
        "bliss": "synergy_bliss", "loewe": "synergy_loewe",
        "hsa": "synergy_hsa", "zip": "synergy_zip",
    }
    df = pd.read_csv(csv_path)

    missing_cols = [c for c in ("drug_row", "drug_col", "cell_line_name") if c not in df.columns]
    if missing_cols:
        raise ValueError(
            f"Expected columns {missing_cols} not found. Real columns are: "
            f"{df.columns.tolist()}. Run inspect_drugcomb_csv() and adjust this "
            f"loader's column-name assumptions -- do not guess further."
        )
    missing_metric_cols = [c for c in synergy_columns.values() if c not in df.columns]
    if missing_metric_cols:
        raise ValueError(
            f"Expected synergy metric columns {missing_metric_cols} not found. "
            f"Real columns are: {df.columns.tolist()}."
        )

    if cell_lines is not None:
        before = len(df)
        df = df[df["cell_line_name"].isin(cell_lines)]
        print(f"Restricted to {len(cell_lines)} requested cell lines: {before} -> {len(df)} rows.")

    has_smiles_row = df["drug_row"].isin(drug_name_to_smiles)
    has_smiles_col = df["drug_col"].isin(drug_name_to_smiles)
    before = len(df)
    all_drugs_before_filter = set(pd.concat([df["drug_row"], df["drug_col"]]))
    df = df[has_smiles_row & has_smiles_col]
    dropped = before - len(df)
    if dropped:
        missing_drugs = sorted(all_drugs_before_filter - set(drug_name_to_smiles))
        print(f"Dropped {dropped} of {before} rows (one or both drugs lack a SMILES "
              f"mapping -- likely biologics or a naming mismatch, not real gaps in "
              f"DrugComb itself). Check whether these were expected: {missing_drugs[:20]}"
              f"{'...' if len(missing_drugs) > 20 else ''}")

    return df.rename(columns={v: f"synergy_{k}" for k, v in synergy_columns.items()})


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary-csv", required=True)
    parser.add_argument("--inspect-only", action="store_true")
    args = parser.parse_args()
    if args.inspect_only:
        inspect_drugcomb_csv(args.summary_csv)
        return
    raise SystemExit(
        "This is a library module for use from your own script once you have a real "
        "drug-name -> SMILES mapping -- see load_drugcomb_summary()'s docstring. "
        "Run with --inspect-only to check column names against your real file first."
    )


if __name__ == "__main__":
    main()
