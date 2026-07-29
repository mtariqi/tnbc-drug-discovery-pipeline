from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import torch
from torch.utils.data import DataLoader

from src.predictive_models.data.dataset import PerturbationArrays, PerturbationDataset
from src.predictive_models.data.splits import make_splits
from src.predictive_models.metrics import evaluate_arrays
from src.predictive_models.models import ConditionalVectorField
from src.predictive_models.train import build_model
from src.predictive_models.utils import load_config, seed_everything, write_json


@torch.no_grad()
def run(config_path: str, checkpoint_path: str) -> dict[str, float]:
    config = load_config(config_path)
    seed_everything(config["seed"])
    arrays = PerturbationArrays.load(config["data"]["path"])
    splits = make_splits(
        arrays.context, arrays.perturbation, config["data"]["split"],
        config["data"]["val_fraction"], config["data"]["test_fraction"], config["seed"]
    )
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = build_model(config, arrays).to(device)
    checkpoint = torch.load(checkpoint_path, map_location=device, weights_only=True)
    model.load_state_dict(checkpoint["model"])
    model.eval()
    loader = DataLoader(PerturbationDataset(arrays, splits.test),
                        batch_size=config["training"]["batch_size"])
    observed, predicted = [], []
    for control, context, perturbation, response in loader:
        control, context, perturbation = control.to(device), context.to(device), perturbation.to(device)
        if isinstance(model, ConditionalVectorField):
            samples = [
                model.sample(control, context, perturbation, config["evaluation"]["flow_steps"])
                for _ in range(config["evaluation"]["samples_per_condition"])
            ]
            estimate = torch.stack(samples).mean(0)
        else:
            estimate = model(control, context, perturbation)
        observed.append(response.numpy())
        predicted.append(estimate.cpu().numpy())
    observed_np, predicted_np = np.concatenate(observed), np.concatenate(predicted)
    metrics = evaluate_arrays(observed_np, predicted_np, config["evaluation"]["de_top_k"])
    output = Path(config["output_dir"])
    write_json(metrics, output / "metrics.json")
    np.savez_compressed(output / "predictions.npz", observed=observed_np, predicted=predicted_np,
                        indices=splits.test)
    print(metrics)
    return metrics


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--checkpoint", required=True)
    args = parser.parse_args()
    run(args.config, args.checkpoint)


if __name__ == "__main__":
    main()

