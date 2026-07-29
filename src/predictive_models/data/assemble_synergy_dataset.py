"""
Assembles the final real-data training file for train_synergy.py, combining:
  1. real_drugcomb_loader.load_drugcomb_summary() -- real synergy records
  2. fetch_smiles.resolve_drug_list() -- drug name -> SMILES mapping
  3. encoders.molecule.featurize_drug_table() -- SMILES -> atom/adjacency arrays

Run as: python -m src.predictive_models.data.assemble_synergy_dataset \\
    --summary-csv <your real DrugComb summary.csv> \\
    --output data/processed/synergy_real.npz

NOT EXECUTED THIS SESSION -- no RDKit, no network access, and no real DrugComb
file were available in the environment this was built in. Every piece was tested
independently (loader logic against a synthetic DataFrame, SMILES featurization
logic reviewed against RDKit's documented API, dataset assembly logic below is
straightforward array-building with no complex control flow) but never run
end-to-end together. Treat your first real run as the actual integration test,
and send me the traceback if anything breaks -- this is exactly the kind of
multi-step real-data glue code most likely to have a small mismatch somewhere.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from src.predictive_models.data.fetch_smiles import resolve_drug_list
from src.predictive_models.data.real_drugcomb_loader import load_drugcomb_summary
from src.predictive_models.encoders.molecule import featurize_drug_table


def assemble(
    summary_csv: str,
    max_atoms: int = 60,
    cell_lines: list[str] | None = None,
    smiles_cache: str = "data/processed/drug_smiles_cache.json",
) -> dict[str, np.ndarray]:
    import pandas as pd

    # Step 1: figure out which drugs actually appear in the summary file, so we
    # only resolve SMILES for drugs we'll actually use (not PubChem's entire catalog).
    raw = pd.read_csv(summary_csv, usecols=lambda c: c in ("drug_row", "drug_col"))
    all_drug_names = sorted(set(raw["drug_row"]) | set(raw["drug_col"]))
    print(f"{len(all_drug_names)} unique drug names found in {summary_csv}.")

    # Step 2: resolve SMILES (PubChem fallback; prefer DrugComb's own metadata file
    # if you have one -- see fetch_smiles.py's docstring).
    drug_name_to_smiles = resolve_drug_list(all_drug_names, cache_path=smiles_cache)

    # Step 3: load the real synergy records, restricted to drugs with a resolved SMILES.
    df = load_drugcomb_summary(summary_csv, drug_name_to_smiles, cell_lines=cell_lines)

    # Step 4: featurize every resolved drug's molecular graph (once per unique drug,
    # not per row).
    used_drugs = sorted(set(df["drug_row"]) | set(df["drug_col"]))
    smiles_subset = {d: drug_name_to_smiles[d] for d in used_drugs}
    featurized = featurize_drug_table(smiles_subset, max_atoms=max_atoms)

    # Drugs that failed featurization (shouldn't happen often at this point, since
    # they already passed the SMILES-resolution step, but a resolved SMILES can
    # still fail RDKit parsing or exceed max_atoms) must also be dropped from df,
    # consistently with how missing-SMILES drugs were dropped in load_drugcomb_summary.
    valid_drugs = set(featurized.keys())
    before = len(df)
    df = df[df["drug_row"].isin(valid_drugs) & df["drug_col"].isin(valid_drugs)]
    if len(df) < before:
        print(f"Dropped {before - len(df)} additional rows where a drug's SMILES "
              f"resolved but failed RDKit featurization or exceeded max_atoms={max_atoms}.")

    drug_list = sorted(valid_drugs)
    drug_to_idx = {d: i for i, d in enumerate(drug_list)}
    cell_line_list = sorted(df["cell_line_name"].unique())
    cell_to_idx = {c: i for i, c in enumerate(cell_line_list)}

    atom_feature_dim = next(iter(featurized.values()))["atom_features"].shape[1]
    drug_atom_features = np.zeros((len(drug_list), max_atoms, atom_feature_dim), dtype=np.float32)
    drug_adjacency = np.zeros((len(drug_list), max_atoms, max_atoms), dtype=np.float32)
    drug_atom_mask = np.zeros((len(drug_list), max_atoms), dtype=np.float32)
    for d, idx in drug_to_idx.items():
        drug_atom_features[idx] = featurized[d]["atom_features"]
        drug_adjacency[idx] = featurized[d]["adjacency"]
        drug_atom_mask[idx] = featurized[d]["atom_mask"]

    data = {
        "drug_row": df["drug_row"].map(drug_to_idx).to_numpy(dtype=np.int64),
        "drug_col": df["drug_col"].map(drug_to_idx).to_numpy(dtype=np.int64),
        "cell_line_id": df["cell_line_name"].map(cell_to_idx).to_numpy(dtype=np.int64),
        "drug_atom_features": drug_atom_features,
        "drug_adjacency": drug_adjacency,
        "drug_atom_mask": drug_atom_mask,
    }
    for metric in ("bliss", "loewe", "hsa", "zip"):
        col = f"synergy_{metric}"
        if col in df.columns:
            data[col] = df[col].to_numpy(dtype=np.float32)
        else:
            print(f"WARNING: {col} not present after loading -- check "
                  f"load_drugcomb_summary's synergy_columns mapping against your real file.")

    print(f"Final assembled dataset: {len(df)} rows, {len(drug_list)} unique drugs, "
          f"{len(cell_line_list)} unique cell lines.")
    return data


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary-csv", required=True)
    parser.add_argument("--output", default="data/processed/synergy_real.npz")
    parser.add_argument("--max-atoms", type=int, default=60)
    args = parser.parse_args()

    data = assemble(args.summary_csv, max_atoms=args.max_atoms)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output, **data)
    print(f"Wrote {output} -- point configs/synergy.yaml's data.path at this file.")


if __name__ == "__main__":
    main()
