"""
MC-Dropout robustness ranking for the synergy GNN.

Runs each candidate drug combination through the trained model N times with
dropout kept ACTIVE at inference (model.train() mode -- safe here since every
normalization layer in this architecture is LayerNorm, not BatchNorm, so there
is no running-statistics pitfall from being in train mode without a real batch).
Produces a DISTRIBUTION of predicted synergy per combination instead of a single
point estimate, then ranks combinations by robustness, not just mean synergy.

IMPORTANT, stated plainly rather than glossed over: this is only meaningful once
the underlying model beats a trivial baseline. Run against the current
DrugComb-v1 checkpoint (mean R2 -0.1374 across all four metrics, per
results/drugcomb_baseline/evaluation_v1.txt) and this will faithfully report
high uncertainty and unstable rankings -- an accurate reflection of that
checkpoint's real state, not a bug in this script, but not yet a useful
'robustness ranking' feature either, since there's no real signal underneath to
be robust about. Point --checkpoint at a retrained model (RTK/NRTK + DepMap
features) once that exists.

Usage:
  python -m src.predictive_models.mc_dropout_robustness \
      --config src/predictive_models/configs/synergy.yaml \
      --checkpoint outputs/synergy/best.pt \
      --n-samples 100 \
      --top-k 3
"""
from __future__ import annotations

import argparse
import json

import numpy as np
import torch

from src.predictive_models.models.synergy_gnn import SynergyGNN
from src.predictive_models.utils import load_config


def mc_dropout_predict(
    model: SynergyGNN,
    drug_atoms: torch.Tensor, drug_adjacency: torch.Tensor, drug_mask: torch.Tensor,
    row_ids: np.ndarray, col_ids: np.ndarray, cell_ids: torch.Tensor,
    n_samples: int,
) -> np.ndarray:
    """Returns predictions of shape [n_samples, n_combinations, n_metrics] in
    NORMALIZED units (caller de-normalizes). Dropout stays active throughout;
    no gradients are computed or needed.
    """
    model.train()  # enables dropout; safe here, see module docstring
    samples = []
    with torch.no_grad():
        for _ in range(n_samples):
            pred = model(
                drug_atoms[row_ids], drug_adjacency[row_ids], drug_mask[row_ids],
                drug_atoms[col_ids], drug_adjacency[col_ids], drug_mask[col_ids],
                cell_ids,
            )
            samples.append(pred.numpy())
    return np.stack(samples, axis=0)


def robustness_scores(
    samples: np.ndarray, metric_mean: np.ndarray, metric_std: np.ndarray,
    metric_idx: int, k: float = 1.0,
) -> dict:
    """samples: [n_samples, n_combinations, n_metrics] in normalized units.
    Returns per-combination mean, std, and the conservative
    (mean - k*std) robustness score, all de-normalized to real units for the
    chosen metric.
    """
    raw = samples[:, :, metric_idx] * metric_std[metric_idx] + metric_mean[metric_idx]
    mean = raw.mean(axis=0)
    std = raw.std(axis=0)
    conservative_score = mean - k * std
    return {"mean": mean, "std": std, "conservative_score": conservative_score}


def rank_stability(samples: np.ndarray, metric_idx: int, top_k: int) -> np.ndarray:
    """For each Monte Carlo sample, rank all combinations by that sample's
    predicted value for this metric; return, per combination, the fraction of
    samples in which it lands in the top_k. This is the rank-stability
    probability -- a combination that's top-k in 95% of samples is a
    genuinely robust recommendation; one that's top-k in 40% of samples is
    not, even if its MEAN rank looks similar.
    """
    n_samples, n_combos, _ = samples.shape
    in_top_k = np.zeros(n_combos, dtype=np.float64)
    for s in range(n_samples):
        ranked = np.argsort(-samples[s, :, metric_idx])  # descending
        in_top_k[ranked[:top_k]] += 1
    return in_top_k / n_samples


def run(config_path: str, checkpoint_path: str, n_samples: int, top_k: int, k: float) -> None:
    config = load_config(config_path)
    arrays = dict(np.load(config["data"]["path"]))
    metrics = config["model"]["synergy_metrics"]

    with open(config["output_dir"] + "/training.json") as f:
        training_info = json.load(f)
    metric_mean = np.array(training_info["metric_mean"])
    metric_std = np.array(training_info["metric_std"])

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
    model.load_state_dict(torch.load(checkpoint_path, map_location="cpu"))

    drug_atoms = torch.from_numpy(arrays["drug_atom_features"]).float()
    drug_adjacency = torch.from_numpy(arrays["drug_adjacency"]).float()
    drug_mask = torch.from_numpy(arrays["drug_atom_mask"]).float()

    # Default candidate set: every unique (drug_row, drug_col, cell_line) triple
    # in the loaded data. Swap this out for your actual focal-patient candidate
    # list once you have real drug-target/RTK-NRTK-derived candidates to rank.
    row_ids = arrays["drug_row"]
    col_ids = arrays["drug_col"]
    cell_ids = torch.from_numpy(arrays["cell_line_id"]).long()

    print(f"Running {n_samples} MC-Dropout samples over {len(row_ids)} candidate rows...")
    samples = mc_dropout_predict(model, drug_atoms, drug_adjacency, drug_mask, row_ids, col_ids, cell_ids, n_samples)

    for metric_idx, metric_name in enumerate(metrics):
        scores = robustness_scores(samples, metric_mean, metric_std, metric_idx, k=k)
        stability = rank_stability(samples, metric_idx, top_k)

        order = np.argsort(-scores["conservative_score"])[:10]
        print(f"\n=== {metric_name}: top 10 by conservative score (mean - {k}*std) ===")
        print(f"{'row':>6}{'col':>6}{'cell':>6}{'mean':>10}{'std':>10}{'cons.score':>12}{'top'+str(top_k)+'_stability':>16}")
        for i in order:
            print(f"{row_ids[i]:>6}{col_ids[i]:>6}{arrays['cell_line_id'][i]:>6}"
                  f"{scores['mean'][i]:>10.3f}{scores['std'][i]:>10.3f}"
                  f"{scores['conservative_score'][i]:>12.3f}{stability[i]:>16.2f}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--checkpoint", default=None, help="Defaults to <output_dir>/best.pt from config")
    parser.add_argument("--n-samples", type=int, default=100)
    parser.add_argument("--top-k", type=int, default=3)
    parser.add_argument("--k", type=float, default=1.0, help="Std-dev penalty weight for conservative score")
    args = parser.parse_args()

    config = load_config(args.config)
    checkpoint = args.checkpoint or (config["output_dir"] + "/best.pt")
    run(args.config, checkpoint, args.n_samples, args.top_k, args.k)
