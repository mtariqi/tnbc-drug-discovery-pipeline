"""
Training entry point: `python -m src.predictive_models.train_synergy --config configs/synergy.yaml`

Uses the dedicated leakage-safe unseen_drug_split for the "unseen_perturbation"-style
evaluation (see data/synergy_splits.py's docstring for why the generic splits.py
can't be reused naively here), and the existing splits.make_splits directly for
unseen_context (cell line), which has no such ambiguity.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import torch
from torch.utils.data import DataLoader, Dataset

from src.predictive_models.data.splits import make_splits
from src.predictive_models.data.synergy_splits import unseen_drug_split
from src.predictive_models.models.synergy_gnn import SynergyGNN
from src.predictive_models.utils import load_config, seed_everything, write_json


class SynergyDataset(Dataset):
    def __init__(self, arrays: dict[str, np.ndarray], indices: np.ndarray, metrics: list[str],
                 metric_mean: np.ndarray, metric_std: np.ndarray):
        self.arrays = arrays
        self.indices = indices
        self.metrics = metrics
        self.metric_mean = torch.from_numpy(metric_mean).float()
        self.metric_std = torch.from_numpy(metric_std).float()
        self.drug_atoms = torch.from_numpy(arrays["drug_atom_features"]).float()
        self.drug_adjacency = torch.from_numpy(arrays["drug_adjacency"]).float()
        self.drug_mask = torch.from_numpy(arrays["drug_atom_mask"]).float()

    def __len__(self) -> int:
        return len(self.indices)

    def __getitem__(self, item: int) -> dict:
        idx = self.indices[item]
        row_id = int(self.arrays["drug_row"][idx])
        col_id = int(self.arrays["drug_col"][idx])
        raw_target = torch.tensor(
            [self.arrays[f"synergy_{m}"][idx] for m in self.metrics], dtype=torch.float32
        )
        # Standardize each metric using TRAINING-set statistics (passed in, computed
        # once in run() from the train split only, never from val/test -- avoids
        # leakage). Without this, metrics with a much larger raw scale (e.g. Loewe,
        # ~6x the variance of the others in this dataset) dominate the combined MSE
        # loss almost completely -- confirmed as the likely cause of a real training
        # run underperforming a trivial mean-baseline (val_mse 103.6 vs baseline 74.3).
        target = (raw_target - self.metric_mean) / self.metric_std
        return {
            "drug_row_atoms": self.drug_atoms[row_id], "drug_row_adjacency": self.drug_adjacency[row_id],
            "drug_row_mask": self.drug_mask[row_id],
            "drug_col_atoms": self.drug_atoms[col_id], "drug_col_adjacency": self.drug_adjacency[col_id],
            "drug_col_mask": self.drug_mask[col_id],
            "cell_line_id": torch.tensor(int(self.arrays["cell_line_id"][idx]), dtype=torch.long),
            "target": target,
        }


def _epoch(model, loader, optimizer, device) -> float:
    training = optimizer is not None
    model.train(training)
    total_loss, n = 0.0, 0
    for batch in loader:
        batch = {k: v.to(device) for k, v in batch.items()}
        with torch.set_grad_enabled(training):
            pred = model(
                batch["drug_row_atoms"], batch["drug_row_adjacency"], batch["drug_row_mask"],
                batch["drug_col_atoms"], batch["drug_col_adjacency"], batch["drug_col_mask"],
                batch["cell_line_id"],
            )
            loss = torch.nn.functional.mse_loss(pred, batch["target"])
            if training:
                optimizer.zero_grad(set_to_none=True)
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
                optimizer.step()
        batch_size = batch["target"].shape[0]
        total_loss += loss.item() * batch_size
        n += batch_size
    return total_loss / n


def run(config_path: str) -> None:
    config = load_config(config_path)
    seed_everything(config["seed"])
    Path(config["output_dir"]).mkdir(parents=True, exist_ok=True)

    arrays = dict(np.load(config["data"]["path"]))
    metrics = config["model"]["synergy_metrics"]

    split_type = config["data"]["split"]
    if split_type == "unseen_drug":
        split = unseen_drug_split(
            arrays["drug_row"], arrays["drug_col"],
            config["data"]["val_fraction"], config["data"]["test_fraction"], config["seed"],
        )
    elif split_type == "unseen_context":
        split = make_splits(
            arrays["cell_line_id"], arrays["drug_row"], "unseen_context",
            config["data"]["val_fraction"], config["data"]["test_fraction"], config["seed"],
        )
    else:
        raise ValueError(f"Unknown split type {split_type!r} -- use 'unseen_drug' or 'unseen_context'.")

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    n_cell_lines = int(arrays["cell_line_id"].max()) + 1

    # Per-metric standardization stats, computed from TRAINING rows only (never
    # val/test, to avoid leakage). See SynergyDataset.__getitem__ for why this
    # matters -- without it, Loewe's much larger raw scale dominated the loss.
    metric_mean = np.array([arrays[f"synergy_{m}"][split.train].mean() for m in metrics], dtype=np.float32)
    metric_std = np.array([arrays[f"synergy_{m}"][split.train].std() for m in metrics], dtype=np.float32)
    metric_std = np.where(metric_std < 1e-6, 1.0, metric_std)  # guard against a degenerate all-constant metric
    print(f"Per-metric training stats -- mean: {dict(zip(metrics, metric_mean.tolist()))}, "
          f"std: {dict(zip(metrics, metric_std.tolist()))}")

    model = SynergyGNN(
        atom_feature_dim=arrays["drug_atom_features"].shape[-1],
        n_cell_lines=n_cell_lines,
        latent_dim=config["model"]["latent_dim"],
        gnn_layers=config["model"]["gnn_layers"],
        cell_embedding_dim=config["model"]["cell_embedding_dim"],
        synergy_metrics=tuple(metrics),
        dropout=config["model"]["dropout"],
    ).to(device)

    def make_loader(indices, shuffle):
        return DataLoader(
            SynergyDataset(arrays, indices, metrics, metric_mean, metric_std),
            batch_size=config["training"]["batch_size"], shuffle=shuffle,
        )

    train_loader = make_loader(split.train, True)
    val_loader = make_loader(split.val, False)

    optimizer = torch.optim.AdamW(
        model.parameters(), lr=config["training"]["learning_rate"], weight_decay=config["training"]["weight_decay"]
    )

    best, stale = float("inf"), 0
    history = []
    for epoch in range(config["training"]["epochs"]):
        train_loss = _epoch(model, train_loader, optimizer, device)
        val_loss = _epoch(model, val_loader, None, device)
        history.append({"epoch": epoch + 1, "train_mse": train_loss, "val_mse": val_loss})
        print(f"epoch {epoch+1}: train_mse={train_loss:.4f} val_mse={val_loss:.4f}")
        if val_loss < best:
            best, stale = val_loss, 0
            torch.save(model.state_dict(), config["output_dir"] + "/best.pt")
        else:
            stale += 1
            if stale >= config["training"]["patience"]:
                break

    write_json(
        {
            "best_val_mse": best, "history": history,
            "metrics": metrics, "metric_mean": metric_mean.tolist(), "metric_std": metric_std.tolist(),
            "split_type": split_type,
        },
        config["output_dir"] + "/training.json",
    )


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    run(args.config)
