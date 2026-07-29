"""
Synthetic data generator shaped like DrugComb's real summary file, for smoke-testing
train_synergy.py without RDKit or a real DrugComb download. Produces fake but
internally-consistent "molecular graphs" (random atom features/adjacency standing
in for RDKit-featurized SMILES) so the model's shapes and training loop can be
verified end to end.

Swap for real_drugcomb_loader.py's output once you have your real DrugComb pull
and RDKit available -- see that module's docstring for the exact real column names
this needs to be verified against.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def generate(
    n_pairs: int = 4000,
    n_drugs: int = 60,
    n_cell_lines: int = 30,
    max_atoms: int = 60,
    atom_feature_dim: int = 21,  # matches encoders/molecule.py's ATOM_FEATURE_DIM
    seed: int = 11,
) -> dict[str, np.ndarray]:
    rng = np.random.default_rng(seed)

    # One synthetic "molecular graph" per drug, reused across every pair it
    # appears in -- mirrors how real drugs are featurized once, not per-row.
    drug_atom_features = np.zeros((n_drugs, max_atoms, atom_feature_dim), dtype=np.float32)
    drug_adjacency = np.zeros((n_drugs, max_atoms, max_atoms), dtype=np.float32)
    drug_atom_mask = np.zeros((n_drugs, max_atoms), dtype=np.float32)
    # each drug's "molecular latent identity" -- a random low-dim vector used only
    # to make synergy labels below depend coherently on which drugs are paired,
    # not to represent real chemistry
    drug_latent = rng.normal(size=(n_drugs, 8)).astype(np.float32)

    for d in range(n_drugs):
        n_atoms = rng.integers(15, max_atoms)
        drug_atom_mask[d, :n_atoms] = 1.0
        drug_atom_features[d, :n_atoms] = rng.normal(size=(n_atoms, atom_feature_dim))
        # random sparse-ish molecular graph
        for i in range(n_atoms):
            n_edges = rng.integers(1, 3)
            partners = rng.choice([j for j in range(n_atoms) if j != i], size=min(n_edges, n_atoms - 1), replace=False)
            drug_adjacency[d, i, partners] = 1.0
            drug_adjacency[d, partners, i] = 1.0

    cell_latent = rng.normal(size=(n_cell_lines, 8)).astype(np.float32)

    drug_row = rng.integers(0, n_drugs, n_pairs)
    drug_col = rng.integers(0, n_drugs, n_pairs)
    cell_line_id = rng.integers(0, n_cell_lines, n_pairs)

    # Synthetic ground truth: synergy metrics are a fixed (but nontrivial, symmetric
    # in drug order) function of the two drugs' latent vectors and the cell line's
    # latent vector, plus noise -- gives the model something real to learn and lets
    # you check that a trained model beats a random/mean baseline.
    metrics = {}
    for name, seed_offset in [("bliss", 1), ("loewe", 2), ("hsa", 3), ("zip", 4)]:
        proj = rng.normal(size=8, scale=0.5).astype(np.float32)
        row_l, col_l, cell_l = drug_latent[drug_row], drug_latent[drug_col], cell_latent[cell_line_id]
        symmetric_term = (row_l + col_l) @ proj + (row_l * col_l) @ proj + cell_l @ proj
        metrics[f"synergy_{name}"] = (symmetric_term + rng.normal(scale=0.3, size=n_pairs)).astype(np.float32)

    return {
        "drug_row": drug_row.astype(np.int64),
        "drug_col": drug_col.astype(np.int64),
        "cell_line_id": cell_line_id.astype(np.int64),
        "drug_atom_features": drug_atom_features,
        "drug_adjacency": drug_adjacency,
        "drug_atom_mask": drug_atom_mask,
        **metrics,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="data/processed/synergy_synthetic.npz")
    parser.add_argument("--seed", type=int, default=11)
    args = parser.parse_args()
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output, **generate(seed=args.seed))
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
