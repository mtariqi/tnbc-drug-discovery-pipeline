"""
Test-set evaluation: trivial baseline vs. XGBoost vs. the trained GNN, all on the
SAME leakage-safe unseen_drug test split -- directly answers the question the
GNN's own negative validation result raised: is 118 unique drugs too few for ANY
model to learn transferable synergy signal, or does the GNN specifically need
more data/tuning that a simpler model wouldn't?

Usage: python -m src.predictive_models.evaluate_synergy --config src/predictive_models/configs/synergy.yaml

Requires: pip install xgboost rdkit (rdkit should already be installed)
"""
from __future__ import annotations

import argparse
import json

import numpy as np
import torch
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

from src.predictive_models.data.synergy_splits import unseen_drug_split
from src.predictive_models.models.synergy_gnn import SynergyGNN
from src.predictive_models.utils import load_config


def morgan_fingerprints(atom_features: np.ndarray, adjacency: np.ndarray, atom_mask: np.ndarray) -> np.ndarray:
    """XGBoost needs fixed-length feature vectors, not molecular graphs. Rather than
    re-deriving Morgan fingerprints from scratch (would need the original SMILES,
    not saved per-drug in the assembled .npz), this uses a graph-level summary
    statistic of the ALREADY-COMPUTED atom features as the feature vector: masked
    mean and max pooling over real (non-padding) atoms. This is a legitimate,
    cheap graph descriptor -- not as expressive as a real circular fingerprint, but
    fair to compare against, since it's derived from the same underlying atom
    featurization the GNN itself uses, not extra chemical information the GNN
    lacks. If this baseline still meaningfully beats the GNN, that is a real
    signal the GNN specifically (not just molecular representation in general)
    needs more data.
    """
    mask = atom_mask[..., None]  # [n_drugs, max_atoms, 1]
    masked = atom_features * mask
    counts = mask.sum(axis=1).clip(min=1)
    mean_pool = masked.sum(axis=1) / counts
    max_pool = np.where(mask.astype(bool), atom_features, -np.inf).max(axis=1)
    max_pool = np.nan_to_num(max_pool, neginf=0.0)
    return np.concatenate([mean_pool, max_pool], axis=1)  # [n_drugs, 2*feature_dim]


def evaluate_split(y_true: np.ndarray, y_pred: np.ndarray, label: str) -> dict:
    return {
        "label": label,
        "r2": r2_score(y_true, y_pred),
        "mae": mean_absolute_error(y_true, y_pred),
        "mse": mean_squared_error(y_true, y_pred),
    }


def run(config_path: str) -> None:
    config = load_config(config_path)
    arrays = dict(np.load(config["data"]["path"]))
    metrics = config["model"]["synergy_metrics"]

    with open(config["output_dir"] + "/training.json") as f:
        training_info = json.load(f)
    metric_mean = np.array(training_info["metric_mean"])
    metric_std = np.array(training_info["metric_std"])

    # Same split function, same seed -> identical partition to what training used,
    # without needing to have saved the indices explicitly (deterministic given the
    # same inputs).
    split = unseen_drug_split(
        arrays["drug_row"], arrays["drug_col"],
        config["data"]["val_fraction"], config["data"]["test_fraction"], config["seed"],
    )
    test_idx = split.test
    print(f"Evaluating on {len(test_idx)} test rows (unseen drugs, never in train or val).")

    raw_targets = np.stack([arrays[f"synergy_{m}"][test_idx] for m in metrics], axis=1)

    results = []

    # --- Trivial baseline: always predict the training mean ---
    trivial_pred = np.tile(metric_mean, (len(test_idx), 1))
    for i, m in enumerate(metrics):
        results.append(evaluate_split(raw_targets[:, i], trivial_pred[:, i], f"trivial_{m}"))

    # --- XGBoost baseline ---
    try:
        from xgboost import XGBRegressor
    except ImportError:
        print("xgboost not installed (pip install xgboost) -- skipping that baseline, "
              "still reporting trivial baseline and GNN below.")
        XGBRegressor = None

    if XGBRegressor is not None:
        fp = morgan_fingerprints(arrays["drug_atom_features"], arrays["drug_adjacency"], arrays["drug_atom_mask"])
        row_features = fp[arrays["drug_row"]]
        col_features = fp[arrays["drug_col"]]
        cell_onehot = np.eye(int(arrays["cell_line_id"].max()) + 1)[arrays["cell_line_id"]]
        X = np.concatenate([row_features, col_features, cell_onehot], axis=1)

        X_train, X_test = X[split.train], X[test_idx]
        for i, m in enumerate(metrics):
            y_train = arrays[f"synergy_{m}"][split.train]
            y_test = arrays[f"synergy_{m}"][test_idx]
            reg = XGBRegressor(n_estimators=200, max_depth=6, learning_rate=0.05, n_jobs=-1)
            reg.fit(X_train, y_train)
            y_pred = reg.predict(X_test)
            results.append(evaluate_split(y_test, y_pred, f"xgboost_{m}"))

    # --- GNN ---
    n_cell_lines = int(arrays["cell_line_id"].max()) + 1
    model = SynergyGNN(
        atom_feature_dim=arrays["drug_atom_features"].shape[-1],
        n_cell_lines=n_cell_lines,
        latent_dim=config["model"]["latent_dim"],
        gnn_layers=config["model"]["gnn_layers"],
        cell_embedding_dim=config["model"]["cell_embedding_dim"],
        synergy_metrics=tuple(metrics),
        dropout=config["model"]["dropout"],
    )
    model.load_state_dict(torch.load(config["output_dir"] + "/best.pt", map_location="cpu"))
    model.eval()

    drug_atoms = torch.from_numpy(arrays["drug_atom_features"]).float()
    drug_adjacency = torch.from_numpy(arrays["drug_adjacency"]).float()
    drug_mask = torch.from_numpy(arrays["drug_atom_mask"]).float()

    row_ids = arrays["drug_row"][test_idx]
    col_ids = arrays["drug_col"][test_idx]
    cell_ids = torch.from_numpy(arrays["cell_line_id"][test_idx]).long()

    with torch.no_grad():
        pred_normalized = model(
            drug_atoms[row_ids], drug_adjacency[row_ids], drug_mask[row_ids],
            drug_atoms[col_ids], drug_adjacency[col_ids], drug_mask[col_ids],
            cell_ids,
        ).numpy()
    pred_raw = pred_normalized * metric_std + metric_mean  # de-normalize back to real units

    for i, m in enumerate(metrics):
        results.append(evaluate_split(raw_targets[:, i], pred_raw[:, i], f"gnn_{m}"))

    print(f"\n{'label':<20}{'R2':>10}{'MAE':>10}{'MSE':>10}")
    for r in results:
        print(f"{r['label']:<20}{r['r2']:>10.4f}{r['mae']:>10.4f}{r['mse']:>10.4f}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    run(args.config)
