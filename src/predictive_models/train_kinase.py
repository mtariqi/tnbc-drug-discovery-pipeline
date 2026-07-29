"""
Training entry point for KinaseCellModel: `python -m src.predictive_models.train_kinase --config configs/kinase_cell.yaml`

Wires the three encoders + fusion model (previously untrained, see src/predictive_models/models/
kinase_cell.py) to real losses and reuses src.predictive_models.data.splits.make_splits for
leakage-safe OOD evaluation -- no new splitting logic was needed since context
(cell line) and perturbation (gene) group arrays are the same shape/semantics the
existing splits module already expects.

Losses, one per model head:
  - response_mean / response_log_variance: heteroscedastic Gaussian NLL rather than
    plain MSE. This uses the log_variance head KinaseCellModel already outputs but
    that nothing previously trained -- and it matters here specifically because
    combinatorial perturbations are expected to be noisier/less certain than single
    perturbations, which a fixed-variance loss (plain MSE) can't represent.
  - redundancy_score: binary cross-entropy against a redundancy label in [0, 1].
    Replace the synthetic proxy label (see data/kinase_synthetic.py) with real
    single-vs-combination-knockout viability comparisons when available -- the loss
    function doesn't change, only what you feed it as `redundancy_label`.
  - pathway_scores: MSE against a pathway-activity target. Swap for whatever your
    real pathway-activity signal's natural loss is (e.g. BCE if you binarize
    "pathway active/inactive" instead of a continuous score).

The three losses are summed with configurable weights (config["training"]["loss_weights"])
rather than hardcoded 1:1:1, since there's no reason to expect the three heads need
equal weight and you'll likely want to tune this once real labels are in place.

Known limitation carried over from the conditional MLP/flow model review: gene and
context embeddings are learned lookup tables sized to cover every id in the dataset,
including ids that only appear in the held-out OOD split. Those specific embedding
rows never receive a gradient update during training, so at evaluation time an
unseen gene/cell-line's embedding is whatever it was randomly initialized to --
carrying no real information. This doesn't break the training loop, but it caps how
much `unseen_perturbation`/`unseen_context` results can mean without a side-information
encoder (e.g. feeding the kinase's real node features from KinaseNetworkEncoder INTO
the perturbation embedding for genes in the panel, rather than relying purely on a
learned id lookup) -- worth prioritizing before reporting OOD numbers as meaningful.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import torch
from torch.utils.data import DataLoader, Dataset

from src.predictive_models.data.splits import make_splits
from src.predictive_models.encoders import CellStateEncoder, KinaseNetworkEncoder, PerturbationEncoder
from src.predictive_models.models.kinase_cell import KinaseCellModel
from src.predictive_models.utils import load_config, seed_everything, write_json


class KinaseDataset(Dataset):
    """Wraps the flat-array format written by data/kinase_synthetic.py (or a real
    loader producing the same keys) into per-sample dicts the training loop consumes.

    `active_heads` controls which optional labels are read: "redundancy_label" and
    "pathway_scores" are only fetched (and don't need to exist in `arrays` at all)
    if the corresponding head is active -- see compute_loss's docstring for why
    this matters on real data, where neither label currently exists.
    """

    def __init__(
        self, arrays: dict[str, np.ndarray], indices: np.ndarray, modality_names: list[str],
        active_heads: set[str],
    ):
        self.arrays = arrays
        self.indices = indices
        self.modality_names = modality_names
        self.active_heads = active_heads

    def __len__(self) -> int:
        return len(self.indices)

    def __getitem__(self, item: int) -> dict:
        idx = self.indices[item]
        sample = {
            "gene_ids": torch.from_numpy(self.arrays["gene_ids"][idx]),
            "direction_ids": torch.tensor(self.arrays["direction_ids"][idx], dtype=torch.long),
            "dose": torch.tensor(self.arrays["dose"][idx], dtype=torch.float32),
        }
        if "response" in self.active_heads:
            sample["response"] = torch.from_numpy(np.atleast_1d(self.arrays["response"][idx]))
        if "redundancy" in self.active_heads:
            sample["redundancy_label"] = torch.tensor(self.arrays["redundancy_label"][idx], dtype=torch.float32)
        if "pathway" in self.active_heads:
            sample["pathway_scores"] = torch.from_numpy(self.arrays["pathway_scores"][idx])
        for name in self.modality_names:
            sample[f"modality__{name}"] = torch.from_numpy(self.arrays[f"modality__{name}"][idx])
            sample[f"modality_mask__{name}"] = torch.tensor(
                self.arrays[f"modality_mask__{name}"][idx], dtype=torch.float32
            )
        return sample


def collate(batch: list[dict], modality_names: list[str]) -> dict:
    stacked = {key: torch.stack([sample[key] for sample in batch]) for key in batch[0]}
    cell_modalities = {name: stacked.pop(f"modality__{name}") for name in modality_names}
    modality_mask = {name: stacked.pop(f"modality_mask__{name}") for name in modality_names}
    stacked["cell_modalities"] = cell_modalities
    stacked["modality_mask"] = modality_mask
    return stacked


def build_model(config: dict, arrays: dict[str, np.ndarray], modality_names: list[str]) -> KinaseCellModel:
    latent_dim = config["model"]["latent_dim"]
    dropout = config["model"]["dropout"]

    cell_encoder = CellStateEncoder(
        modality_dims={name: arrays[f"modality__{name}"].shape[1] for name in modality_names},
        latent_dim=latent_dim,
        dropout=dropout,
    )
    n_genes = int(arrays["gene_ids"].max()) + 1  # +1 since 0 is the padding index
    perturbation_encoder = PerturbationEncoder(
        n_genes=n_genes,
        n_directions=len(config["perturbation"]["directions"]),
        latent_dim=latent_dim,
        dropout=dropout,
    )
    adjacency = torch.from_numpy(arrays["kinase_adjacency"]).float()
    network_encoder = KinaseNetworkEncoder(
        node_feature_dim=arrays["kinase_node_features"].shape[1],
        adjacency=adjacency,
        latent_dim=latent_dim,
        layers=config["model"]["network_layers"],
        dropout=dropout,
    )
    return KinaseCellModel(
        cell_encoder=cell_encoder,
        perturbation_encoder=perturbation_encoder,
        network_encoder=network_encoder,
        latent_dim=latent_dim,
        response_dim=config["model"]["response_dim"],
        pathway_dim=config["model"]["pathway_dim"],
        dropout=dropout,
    )


def gaussian_nll(target: torch.Tensor, mean: torch.Tensor, log_variance: torch.Tensor) -> torch.Tensor:
    """Heteroscedastic Gaussian negative log-likelihood, elementwise-mean reduced.
    log_variance is already clamped to [-8, 6] inside KinaseCellModel to keep this
    numerically stable (prevents variance collapsing to ~0 and exploding the loss).
    """
    precision = torch.exp(-log_variance)
    return torch.mean(0.5 * (precision * (target - mean) ** 2 + log_variance))


def compute_loss(
    outputs: dict, batch: dict, loss_weights: dict[str, float], active_heads: set[str]
) -> tuple[torch.Tensor, dict[str, float]]:
    """Only heads in `active_heads` contribute to the loss and are required in
    `batch`. This matters for real data: as of this project's own docs/limitations.md,
    no real label exists yet for redundancy or pathway activity (no combination-
    synergy data, no pathway-activity ground truth) -- training those heads would
    mean inventing a proxy label, which this project's own stated discipline says
    not to do. Set their loss_weights to 0 in the config and they're skipped here
    entirely (including not requiring the corresponding array in the dataset),
    rather than silently trained against a placeholder.
    """
    total = torch.zeros((), device=outputs["response_mean"].device)
    components: dict[str, float] = {}

    if "response" in active_heads:
        response_loss = gaussian_nll(batch["response"], outputs["response_mean"], outputs["response_log_variance"])
        total = total + loss_weights["response"] * response_loss
        components["response_nll"] = response_loss.item()

    if "redundancy" in active_heads:
        redundancy_loss = torch.nn.functional.binary_cross_entropy(
            outputs["redundancy_score"].clamp(1e-6, 1 - 1e-6), batch["redundancy_label"]
        )
        total = total + loss_weights["redundancy"] * redundancy_loss
        components["redundancy_bce"] = redundancy_loss.item()

    if "pathway" in active_heads:
        pathway_loss = torch.nn.functional.mse_loss(outputs["pathway_scores"], batch["pathway_scores"])
        total = total + loss_weights["pathway"] * pathway_loss
        components["pathway_mse"] = pathway_loss.item()

    components["total"] = total.item()
    return total, components


def build_target_mask(gene_ids: torch.Tensor, n_kinases: int) -> torch.Tensor:
    """Flag, per sample, which kinase-graph node(s) correspond to the perturbed
    gene(s), so KinaseNetworkEncoder's final attention pool reads out the
    (multi-hop, network-contextualized) representation of the SPECIFIC kinase(s)
    being perturbed -- rather than an identical whole-graph summary for every
    sample regardless of what was perturbed.

    ASSUMPTION (synthetic data only): gene id `g` (1-indexed, 0 = padding) maps
    directly to kinase-panel node index `g - 1`. This holds in
    data/kinase_synthetic.py by construction (n_genes == n_kinases == 90) but will
    NOT generally hold on real data, where perturbed genes can fall outside the
    90-kinase panel entirely (e.g. TP53, a drug target not in the RTK/NRTK panel).
    Replace this with a real gene-symbol -> kinase-panel-index lookup (None/all-
    False for genes outside the panel) once wiring in real perturbation metadata.
    Rows with no valid in-panel target fall back to uniform attention over the
    whole network (all-True), which is a safe default, not a bug: it means "no
    specific kinase to highlight, use the whole-network summary" -- exactly the
    right fallback for panel-external perturbations like TP53.
    """
    batch_size = gene_ids.shape[0]
    mask = torch.zeros(batch_size, n_kinases, dtype=torch.bool, device=gene_ids.device)
    for row in range(batch_size):
        node_indices = [int(g) - 1 for g in gene_ids[row].tolist() if 0 < int(g) <= n_kinases]
        for idx in node_indices:
            mask[row, idx] = True
    no_target = ~mask.any(dim=1)
    mask[no_target] = True  # fall back to attend over the whole graph
    return mask


def _epoch(model, loader, node_features, optimizer, device, loss_weights, active_heads) -> dict[str, float]:
    training = optimizer is not None
    model.train(training)
    totals = {"total": 0.0}
    for head in active_heads:
        totals[{"response": "response_nll", "redundancy": "redundancy_bce", "pathway": "pathway_mse"}[head]] = 0.0
    n = 0
    n_kinases = node_features.shape[0]
    for batch in loader:
        cell_modalities = {k: v.to(device) for k, v in batch["cell_modalities"].items()}
        modality_mask = {k: v.to(device) for k, v in batch["modality_mask"].items()}
        gene_ids = batch["gene_ids"].to(device)
        direction_ids = batch["direction_ids"].to(device)
        dose = batch["dose"].to(device)
        batch_size = gene_ids.shape[0]

        loss_batch = {}
        if "response" in active_heads:
            loss_batch["response"] = batch["response"].to(device)
        if "redundancy" in active_heads:
            loss_batch["redundancy_label"] = batch["redundancy_label"].to(device)
        if "pathway" in active_heads:
            loss_batch["pathway_scores"] = batch["pathway_scores"].to(device)

        # Repeat the (fixed, dataset-wide) kinase graph across the batch so each
        # sample gets its OWN attention pool over it -- see build_target_mask's
        # docstring for why a shared, unbatched graph would silently make the
        # network encoder's output identical regardless of what was perturbed.
        batched_node_features = node_features.to(device).unsqueeze(0).expand(batch_size, -1, -1)
        target_mask = build_target_mask(gene_ids, n_kinases)

        with torch.set_grad_enabled(training):
            outputs = model(
                cell_modalities=cell_modalities,
                gene_ids=gene_ids,
                direction_ids=direction_ids,
                dose=dose,
                kinase_node_features=batched_node_features,
                modality_mask=modality_mask,
                target_mask=target_mask,
            )
            loss, components = compute_loss(outputs, loss_batch, loss_weights, active_heads)
            if training:
                optimizer.zero_grad(set_to_none=True)
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
                optimizer.step()

        for key in totals:
            totals[key] += components[key] * batch_size
        n += batch_size

    return {key: value / n for key, value in totals.items()}


def run(config_path: str) -> Path:
    config = load_config(config_path)
    seed_everything(config["seed"])

    data_path = config["data"]["path"]
    arrays = dict(np.load(data_path))
    modality_names = sorted(config["cell_state"]["modalities"])

    loss_weights = config["training"]["loss_weights"]
    active_heads = {head for head, weight in loss_weights.items() if weight > 0}
    if not active_heads:
        raise ValueError(
            "All loss_weights are 0 -- nothing would be trained. Set at least "
            "training.loss_weights.response > 0."
        )
    print(f"Active heads this run: {sorted(active_heads)} "
          f"(heads with weight=0 are skipped entirely -- no label required for them)")

    if "response" in active_heads:
        actual_response_dim = np.atleast_2d(arrays["response"].reshape(len(arrays["response"]), -1)).shape[1]
        configured = config["model"]["response_dim"]
        if actual_response_dim != configured:
            raise ValueError(
                f"config model.response_dim={configured} doesn't match the actual "
                f"response array's last dimension ({actual_response_dim}). For real "
                f"DepMap dependency_prob data this should be 1, not the synthetic "
                f"generator's response_dim=64 default -- check data.path is pointing "
                f"at what you think it is."
            )

    splits = make_splits(
        arrays["context"], arrays["perturbation"], config["data"]["split"],
        config["data"]["val_fraction"], config["data"]["test_fraction"], config["seed"],
    )

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = build_model(config, arrays, modality_names).to(device)
    node_features = torch.from_numpy(arrays["kinase_node_features"]).float()

    def make_loader(indices: np.ndarray, shuffle: bool) -> DataLoader:
        dataset = KinaseDataset(arrays, indices, modality_names, active_heads)
        return DataLoader(
            dataset, batch_size=config["training"]["batch_size"], shuffle=shuffle,
            collate_fn=lambda batch: collate(batch, modality_names),
        )

    train_loader = make_loader(splits.train, shuffle=True)
    val_loader = make_loader(splits.val, shuffle=False)

    optimizer = torch.optim.AdamW(
        model.parameters(), lr=config["training"]["learning_rate"],
        weight_decay=config["training"]["weight_decay"],
    )
    loss_weights = config["training"]["loss_weights"]

    output = Path(config["output_dir"])
    output.mkdir(parents=True, exist_ok=True)
    checkpoint = output / "best.pt"
    best, stale, history = float("inf"), 0, []

    for epoch in range(config["training"]["epochs"]):
        train_metrics = _epoch(model, train_loader, node_features, optimizer, device, loss_weights, active_heads)
        val_metrics = _epoch(model, val_loader, node_features, None, device, loss_weights, active_heads)
        history.append({"epoch": epoch + 1, "train": train_metrics, "val": val_metrics})
        if val_metrics["total"] < best:
            best, stale = val_metrics["total"], 0
            torch.save({"model": model.state_dict(), "config": config}, checkpoint)
        else:
            stale += 1
            if stale >= config["training"]["patience"]:
                break

    write_json(
        {
            "best_val_loss": best, "history": history,
            "split_sizes": {"train": len(splits.train), "val": len(splits.val), "test": len(splits.test)},
        },
        output / "training.json",
    )
    print(f"Best checkpoint: {checkpoint}")
    return checkpoint


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    run(args.config)


if __name__ == "__main__":
    main()
