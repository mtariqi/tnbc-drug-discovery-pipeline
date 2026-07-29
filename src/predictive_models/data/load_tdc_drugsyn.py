"""
Fast path to real-data training: PyTDC's DrugSyn dataset already ships pre-resolved
SMILES per drug (queried from PubChem by the TDC maintainers, not by us) --
skipping fetch_smiles.py's PubChem-lookup step entirely. Confirmed real, current
info (checked this session): `from tdc.multi_pred import DrugSyn; data =
DrugSyn(name='DrugComb')` gives a NCI-60-derived subset (129 drugs, 59 cell
lines, ~297,098 drug-pair-cell-line rows) or `name='OncoPolyPharmacology'` gives
a smaller Merck-derived subset (38 drugs, 39 cell lines, ~23,052 rows) --
useful for a fast first real run before scaling to the full raw DrugComb portal
dump (real_drugcomb_loader.py / assemble_synergy_dataset.py), which has ~1.4M
records but requires the PubChem name-matching step this path avoids.

TDC's column names are CONFIRMED against a real live TDC install this session
(not just inferred): Drug1_ID, Drug2_ID, Cell_Line_ID, CSS, Synergy_ZIP,
Synergy_Bliss, Synergy_Loewe, Synergy_HSA, Drug1, Drug2, CellLine -- where
Drug1/Drug2 are real SMILES strings and Drug1_ID/Drug2_ID are chemical names.
All four synergy metrics are real, independent columns (better than initially
assumed -- no need to replicate a single label across all four slots).

STILL APPLIES: the biologics limitation. TDC's SMILES come from PubChem the same
way our own fetch_smiles.py would have queried it, so trastuzumab and other
biologics will simply be absent from TDC's drug list, not present-but-wrong.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def inspect_tdc_drugsyn(name: str = "DrugComb") -> None:
    """Print the real column names and a few rows before trusting anything
    downstream -- same inspect-first discipline as the rest of this project."""
    from tdc.multi_pred import DrugSyn

    data = DrugSyn(name=name)
    df = data.get_data(format="df")
    print(f"=== TDC DrugSyn(name={name!r}) ===")
    print(f"Columns: {df.columns.tolist()}")
    print(f"Shape: {df.shape}")
    print(df.head())


def load_tdc_drugsyn(
    name: str = "DrugComb",
    max_atoms: int = 60,
    column_map: dict[str, str] | None = None,
) -> dict[str, np.ndarray]:
    """Load a TDC DrugSyn dataset and assemble it into the exact array contract
    train_synergy.py expects. `column_map` lets you override the guessed column
    names if inspect_tdc_drugsyn() shows something different -- default guesses
    follow TDC's documented multi_pred.DrugSyn convention.
    """
    from tdc.multi_pred import DrugSyn

    from src.predictive_models.encoders.molecule import featurize_drug_table

    column_map = column_map or {
        "drug1_id": "Drug1_ID", "drug1_smiles": "Drug1",
        "drug2_id": "Drug2_ID", "drug2_smiles": "Drug2",
        "cell_line": "Cell_Line_ID",
        "bliss": "Synergy_Bliss", "loewe": "Synergy_Loewe",
        "hsa": "Synergy_HSA", "zip": "Synergy_ZIP",
    }

    data = DrugSyn(name=name)
    df = data.get_data(format="df")

    missing = [c for c in column_map.values() if c not in df.columns]
    if missing:
        raise ValueError(
            f"Expected columns {missing} not found. Real columns are: "
            f"{df.columns.tolist()}. Run inspect_tdc_drugsyn() and pass a "
            f"corrected column_map -- do not guess further."
        )

    missing_metric_cols = [c for c in (column_map["bliss"], column_map["loewe"], column_map["hsa"], column_map["zip"]) if c not in df.columns]
    if missing_metric_cols:
        raise ValueError(f"Expected synergy metric columns {missing_metric_cols} not found in {df.columns.tolist()}.")

    before_nan_drop = len(df)
    df = df.dropna(subset=[column_map["bliss"], column_map["loewe"], column_map["hsa"], column_map["zip"]])
    if len(df) < before_nan_drop:
        print(f"Dropped {before_nan_drop - len(df)} rows with a missing (NaN) synergy "
              f"value in at least one of the four metrics -- real biological data "
              f"commonly has gaps here; left in, these would silently corrupt the "
              f"training loss rather than error out.")

    # Build a unique drug -> SMILES table from both drug columns combined
    drug_smiles: dict[str, str] = {}
    for _, row in df.iterrows():
        drug_smiles[row[column_map["drug1_id"]]] = row[column_map["drug1_smiles"]]
        drug_smiles[row[column_map["drug2_id"]]] = row[column_map["drug2_smiles"]]

    print(f"{len(drug_smiles)} unique drugs, {df[column_map['cell_line']].nunique()} "
          f"unique cell lines, {len(df)} pair-rows in TDC's {name!r} dataset.")

    featurized = featurize_drug_table(drug_smiles, max_atoms=max_atoms)
    valid_drugs = set(featurized.keys())
    before = len(df)
    df = df[df[column_map["drug1_id"]].isin(valid_drugs) & df[column_map["drug2_id"]].isin(valid_drugs)]
    if len(df) < before:
        print(f"Dropped {before - len(df)} rows where a drug's SMILES failed RDKit "
              f"featurization or exceeded max_atoms={max_atoms}.")

    drug_list = sorted(valid_drugs)
    drug_to_idx = {d: i for i, d in enumerate(drug_list)}
    cell_line_list = sorted(df[column_map["cell_line"]].unique())
    cell_to_idx = {c: i for i, c in enumerate(cell_line_list)}

    atom_feature_dim = next(iter(featurized.values()))["atom_features"].shape[1]
    drug_atom_features = np.zeros((len(drug_list), max_atoms, atom_feature_dim), dtype=np.float32)
    drug_adjacency = np.zeros((len(drug_list), max_atoms, max_atoms), dtype=np.float32)
    drug_atom_mask = np.zeros((len(drug_list), max_atoms), dtype=np.float32)
    for d, idx in drug_to_idx.items():
        drug_atom_features[idx] = featurized[d]["atom_features"]
        drug_adjacency[idx] = featurized[d]["adjacency"]
        drug_atom_mask[idx] = featurized[d]["atom_mask"]

    # CONFIRMED against a real TDC DrugComb pull this session: unlike our original
    # guess, TDC provides all four synergy metrics as real, independent columns
    # (not one label needing replication) -- genuinely better than assumed.
    data_out = {
        "drug_row": df[column_map["drug1_id"]].map(drug_to_idx).to_numpy(dtype=np.int64),
        "drug_col": df[column_map["drug2_id"]].map(drug_to_idx).to_numpy(dtype=np.int64),
        "cell_line_id": df[column_map["cell_line"]].map(cell_to_idx).to_numpy(dtype=np.int64),
        "drug_atom_features": drug_atom_features,
        "drug_adjacency": drug_adjacency,
        "drug_atom_mask": drug_atom_mask,
        "synergy_bliss": df[column_map["bliss"]].to_numpy(dtype=np.float32),
        "synergy_loewe": df[column_map["loewe"]].to_numpy(dtype=np.float32),
        "synergy_hsa": df[column_map["hsa"]].to_numpy(dtype=np.float32),
        "synergy_zip": df[column_map["zip"]].to_numpy(dtype=np.float32),
    }
    print(f"Final assembled dataset: {len(df)} rows, {len(drug_list)} unique drugs, "
          f"{len(cell_line_list)} unique cell lines. All four synergy_* fields are "
          f"real, independent metrics from TDC's DrugComb data.")
    return data_out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tdc-name", default="DrugComb", choices=["DrugComb", "OncoPolyPharmacology"])
    parser.add_argument("--output", default="data/processed/synergy_tdc.npz")
    parser.add_argument("--max-atoms", type=int, default=60)
    parser.add_argument("--inspect-only", action="store_true")
    args = parser.parse_args()

    if args.inspect_only:
        inspect_tdc_drugsyn(args.tdc_name)
        return

    data = load_tdc_drugsyn(args.tdc_name, max_atoms=args.max_atoms)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output, **data)
    print(f"Wrote {output} -- point configs/synergy.yaml's data.path at this file.")


if __name__ == "__main__":
    main()
