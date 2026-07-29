from __future__ import annotations

import argparse
from pathlib import Path

import torch
from torch.utils.data import DataLoader

from src.predictive_models.data.dataset import PerturbationArrays, PerturbationDataset
from src.predictive_models.data.splits import make_splits
from src.predictive_models.models import ConditionalMLP, ConditionalVectorField
from src.predictive_models.utils import load_config, seed_everything, write_json


def build_model(config, arrays):
    kwargs = dict(
        n_genes=arrays.response.shape[1],
        n_contexts=int(arrays.context.max()) + 1,
        n_perturbations=int(arrays.perturbation.max()) + 1,
        hidden_dim=config["model"]["hidden_dim"],
        context_dim=config["model"]["context_dim"],
        perturbation_dim=config["model"]["perturbation_dim"],
    )
    if config["model"]["type"] == "mlp":
        return ConditionalMLP(**kwargs)
    if config["model"]["type"] == "flow":
        return ConditionalVectorField(**kwargs, time_dim=config["model"]["time_dim"])
    raise ValueError(f"Unknown model type: {config['model']['type']}")


def _epoch(model, loader, optimizer, device):
    training = optimizer is not None
    model.train(training)
    total = 0.0
    for control, context, perturbation, response in loader:
        control, context = control.to(device), context.to(device)
        perturbation, response = perturbation.to(device), response.to(device)
        with torch.set_grad_enabled(training):
            if isinstance(model, ConditionalVectorField):
                loss = model.loss(response, control, context, perturbation)
            else:
                loss = torch.nn.functional.mse_loss(
                    model(control, context, perturbation), response
                )
            if training:
                optimizer.zero_grad(set_to_none=True)
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
                optimizer.step()
        total += loss.item() * len(control)
    return total / len(loader.dataset)


def run(config_path: str) -> Path:
    config = load_config(config_path)
    seed_everything(config["seed"])
    arrays = PerturbationArrays.load(config["data"]["path"])
    splits = make_splits(
        arrays.context, arrays.perturbation, config["data"]["split"],
        config["data"]["val_fraction"], config["data"]["test_fraction"], config["seed"]
    )
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = build_model(config, arrays).to(device)
    loaders = {
        name: DataLoader(PerturbationDataset(arrays, getattr(splits, name)),
                         batch_size=config["training"]["batch_size"],
                         shuffle=name == "train")
        for name in ("train", "val")
    }
    optimizer = torch.optim.AdamW(model.parameters(), lr=config["training"]["learning_rate"],
                                  weight_decay=config["training"]["weight_decay"])
    output = Path(config["output_dir"])
    output.mkdir(parents=True, exist_ok=True)
    checkpoint = output / "best.pt"
    best, stale, history = float("inf"), 0, []
    for epoch in range(config["training"]["epochs"]):
        train_loss = _epoch(model, loaders["train"], optimizer, device)
        val_loss = _epoch(model, loaders["val"], None, device)
        history.append({"epoch": epoch + 1, "train_loss": train_loss, "val_loss": val_loss})
        if val_loss < best:
            best, stale = val_loss, 0
            torch.save({"model": model.state_dict(), "config": config}, checkpoint)
        else:
            stale += 1
            if stale >= config["training"]["patience"]:
                break
    write_json({"best_val_loss": best, "history": history, "split_sizes": {
        "train": len(splits.train), "val": len(splits.val), "test": len(splits.test)
    }}, output / "training.json")
    print(f"Best checkpoint: {checkpoint}")
    return checkpoint


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    run(args.config)


if __name__ == "__main__":
    main()

