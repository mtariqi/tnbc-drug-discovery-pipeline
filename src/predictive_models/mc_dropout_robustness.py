"""
MC-Dropout robustness ranking for the synergy GNN -- Bayesian approximation of
predictive uncertainty via Monte Carlo Dropout (Gal & Ghahramani, 2016), not a
classical pharmacological Monte Carlo simulation. Framing matters: this
quantifies the trained model's own uncertainty about its predictions, not
variability from real-world dosing/PK/PD randomness.

IMPORTANT, stated plainly rather than glossed over: this is only meaningful
once the underlying model beats a trivial baseline. Run against the current
DrugComb-v1 checkpoint (mean R2 -0.1374 across all four metrics, per
results/drugcomb_baseline/evaluation_v1.txt) and this will faithfully report
high uncertainty and unstable rankings -- an accurate reflection of that
checkpoint's real state, not a bug, but not yet a useful robustness-ranking
feature either, since there's no real signal underneath to be robust about.
Point --checkpoint at a retrained model (RTK/NRTK + DepMap features) once
that exists.

Usage:
  python -m src.predictive_models.mc_dropout_robustness \
      --config src/predictive_models/configs/synergy.yaml \
      --checkpoint outputs/synergy/best.pt \
      --n-samples 100 --top-k 3 --max-candidates 2000
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from src.predictive_models.models.synergy_gnn import SynergyGNN
from src.predictive_models.utils import load_config


class RunningStats:
    """Numerically stable online mean/std computation (Welford's algorithm) --
    verified against numpy's direct mean/std on synthetic data before use.
    Streams the computation one MC sample at a time rather than holding the
    full [n_samples, n_combos, n_metrics] array in memory simultaneously."""

    def __init__(self, shape):
        self.n = 0
        self.mean = np.zeros(shape, dtype=np.float64)
        self.M2 = np.zeros(shape, dtype=np.float64)

    def update(self, x: np.ndarray):
        self.n += 1
        delta = x - self.mean
        self.mean += delta / self.n
        delta2 = x - self.mean
        self.M2 += delta * delta2

    @property
    def variance(self):
        if self.n < 2:
            return np.zeros_like(self.mean)
        return self.M2 / (self.n - 1)

    @property
    def std(self):
        return np.sqrt(self.variance)


def enable_mc_dropout(model):
    """Targets ONLY nn.Dropout submodules, leaving everything else (e.g.
    LayerNorm) in eval mode -- more robust than a blanket model.train() if the
    architecture ever gains a BatchNorm-style layer with running statistics."""
    for m in model.modules():
        if isinstance(m, torch.nn.Dropout):
            m.train()


def predict_single_pass(model, drug_atoms, drug_adj, drug_mask, row_ids, col_ids, cell_ids, batch_size, device):
    """Batched forward pass -- bounds peak memory regardless of dataset size.
    The naive unbatched version was found to require ~9GB just for the input
    tensors on a ~230k-row dataset, before any GNN intermediate activations."""
    outputs = []
    with torch.no_grad():
        for start in range(0, len(row_ids), batch_size):
            end = min(start + batch_size, len(row_ids))
            r, c = row_ids[start:end], col_ids[start:end]
            pred = model(
                drug_atoms[r], drug_adj[r], drug_mask[r],
                drug_atoms[c], drug_adj[c], drug_mask[c],
                cell_ids[start:end],
            )
            outputs.append(pred.cpu().numpy())
    return np.concatenate(outputs, axis=0)


def mc_dropout_ranking(
    model, metric_mean, metric_std, drug_atoms, drug_adj, drug_mask,
    row_ids, col_ids, cell_ids, metrics,
    n_samples=100, top_k=3, k_penalty=1.0, batch_size=1024, device="cpu",
):
    n_combos = len(row_ids)
    n_metrics = len(metrics)
    stats = RunningStats((n_combos, n_metrics))
    topk_counts = np.zeros((n_combos, n_metrics))

    model.eval()
    enable_mc_dropout(model)

    print(f"Running MC-Dropout ({n_samples} passes, {n_combos} combinations)")

    for _ in range(n_samples):
        pred = predict_single_pass(model, drug_atoms, drug_adj, drug_mask, row_ids, col_ids, cell_ids, batch_size, device)
        pred_real = pred * metric_std + metric_mean
        stats.update(pred_real)
        for metric_idx in range(n_metrics):
            ranking = np.argsort(-pred_real[:, metric_idx])
            topk_counts[ranking[:top_k], metric_idx] += 1

    topk_prob = topk_counts / n_samples

    results = {}
    for metric_idx, metric_name in enumerate(metrics):
        mean = stats.mean[:, metric_idx]
        std = stats.std[:, metric_idx]
        conservative_score = mean - k_penalty * std
        ci_low = mean - 1.96 * std
        ci_high = mean + 1.96 * std

        df = pd.DataFrame({
            "drug_row": row_ids, "drug_col": col_ids, "cell_line": cell_ids.cpu().numpy(),
            "mean_synergy": mean, "std": std, "ci_low": ci_low, "ci_high": ci_high,
            "conservative_score": conservative_score, "topk_probability": topk_prob[:, metric_idx],
        })
        df = df.sort_values("conservative_score", ascending=False)
        results[metric_name] = df

    return results


def build_model(config, arrays):
    n_cell_lines = int(arrays["cell_line_id"].max()) + 1
    return SynergyGNN(
        atom_feature_dim=arrays["drug_atom_features"].shape[-1],
        n_cell_lines=n_cell_lines,
        latent_dim=config["model"]["latent_dim"],
        gnn_layers=config["model"]["gnn_layers"],
        cell_embedding_dim=config["model"]["cell_embedding_dim"],
        synergy_metrics=tuple(config["model"]["synergy_metrics"]),
        dropout=config["model"]["dropout"],
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--checkpoint", required=True)
    parser.add_argument("--n-samples", type=int, default=100)
    parser.add_argument("--top-k", type=int, default=3)
    parser.add_argument("--k", type=float, default=1.0)
    parser.add_argument("--batch-size", type=int, default=1024)
    parser.add_argument("--output", default="mc_dropout_results")
    parser.add_argument("--max-candidates", type=int, default=2000,
                         help="Randomly subsample to this many candidate rows for a fast first "
                              "run (230k rows x 100 MC samples otherwise = ~23M forward-pass-"
                              "equivalents, batching prevents a crash but not a long runtime). "
                              "Pass a larger value, or omit for the full dataset.")
    args = parser.parse_args()

    config = load_config(args.config)
    arrays = dict(np.load(config["data"]["path"]))

    with open(Path(config["output_dir"]) / "training.json") as f:
        train_info = json.load(f)
    metric_mean = np.array(train_info["metric_mean"])
    metric_std = np.array(train_info["metric_std"])

    model = build_model(config, arrays)
    model.load_state_dict(torch.load(args.checkpoint, map_location="cpu"))

    drug_atoms = torch.tensor(arrays["drug_atom_features"]).float()
    drug_adj = torch.tensor(arrays["drug_adjacency"]).float()
    drug_mask = torch.tensor(arrays["drug_atom_mask"]).float()

    row_ids = arrays["drug_row"]
    col_ids = arrays["drug_col"]
    cell_line_ids = arrays["cell_line_id"]

    if args.max_candidates is not None and args.max_candidates < len(row_ids):
        rng = np.random.default_rng(0)
        keep = rng.choice(len(row_ids), size=args.max_candidates, replace=False)
        row_ids, col_ids, cell_line_ids = row_ids[keep], col_ids[keep], cell_line_ids[keep]
        print(f"Subsampled to {args.max_candidates} of {len(arrays['drug_row'])} candidate rows "
              f"(--max-candidates 0 or a larger value for the full dataset).")

    cell_ids = torch.tensor(cell_line_ids).long()

    results = mc_dropout_ranking(
        model=model, metric_mean=metric_mean, metric_std=metric_std,
        drug_atoms=drug_atoms, drug_adj=drug_adj, drug_mask=drug_mask,
        row_ids=row_ids, col_ids=col_ids, cell_ids=cell_ids,
        metrics=config["model"]["synergy_metrics"],
        n_samples=args.n_samples, top_k=args.top_k, k_penalty=args.k, batch_size=args.batch_size,
    )

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    for metric, df in results.items():
        file_path = output_dir / f"{metric}_robustness.csv"
        df.to_csv(file_path, index=False)
        print("\n" + "=" * 80)
        print(metric)
        print("=" * 80)
        print(df.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
