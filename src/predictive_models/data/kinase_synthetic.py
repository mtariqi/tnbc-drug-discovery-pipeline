"""
Synthetic data generator shaped like the real KinaseCell pipeline's inputs, so
train_kinase.py can be run and verified end to end without your real STRING/DepMap/
TCGA/DGIdb data. Swap this for a real loader once you're ready to point the trainer
at your actual CTS-derived arrays -- the trainer only depends on the array shapes
and names below, not on this generator.

Design choices, so it's clear what's real signal vs. scaffolding:
  - `context` = cell-line id, `perturbation` = perturbed-gene id, matching the
    (context, perturbation) group arrays that src.predictive_models.data.splits.make_splits
    already expects -- no changes to splits.py were needed to reuse it here.
  - Redundancy labels are synthesized as a function of how many OTHER kinases in the
    network remain "active" (a stand-in for network slack / real DepMap co-essentiality
    signal) -- replace with your real single-vs-combination-knockout viability data
    for the redundancy_head to learn something biologically real rather than this
    proxy.
  - Pathway scores are a random linear projection of the true response, so the
    pathway head has a learnable (if synthetic) target to fit during smoke testing.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def generate_kinase_graph(n_kinases: int = 90, seed: int = 0) -> dict[str, np.ndarray]:
    """Stand-in for load_kinase_graph() output in rtk_nrtk_adapters.py -- same
    shapes (node_features [n_kinases, n_node_features], adjacency [n_kinases,
    n_kinases]) so KinaseNetworkEncoder can be constructed identically either way.
    """
    rng = np.random.default_rng(seed)
    is_rtk = (rng.uniform(size=n_kinases) < 0.4).astype(np.float32)  # ~roughly RTK/NRTK split
    numeric_features = rng.normal(size=(n_kinases, 5)).astype(np.float32)  # centrality, essentiality, survival, druggability, cts
    node_features = np.concatenate([numeric_features, is_rtk[:, None]], axis=1)

    # sparse-ish random graph with a handful of strong edges per node, like a
    # thresholded STRING network would produce
    adjacency = np.zeros((n_kinases, n_kinases), dtype=np.float32)
    for i in range(n_kinases):
        n_edges = rng.integers(2, 6)
        partners = rng.choice([j for j in range(n_kinases) if j != i], size=n_edges, replace=False)
        weights = rng.uniform(0.3, 1.0, size=n_edges).astype(np.float32)
        adjacency[i, partners] = np.maximum(adjacency[i, partners], weights)
        adjacency[partners, i] = np.maximum(adjacency[partners, i], weights)

    return {"node_features": node_features, "adjacency": adjacency, "is_rtk": is_rtk}


def generate(
    n_samples: int = 3000,
    n_cell_lines: int = 25,          # matches your real 25 confirmed TNBC lines
    n_genes: int = 90,               # matches your real 90-kinase RTK/NRTK panel
    n_kinases: int = 90,
    n_directions: int = 4,           # padding, activation, knockout/knockdown, inhibition
    n_modality_dims: dict[str, int] | None = None,
    response_dim: int = 1,           # matches real target: DepMap dependency_prob is a scalar per row
    pathway_dim: int = 12,
    max_combo_size: int = 2,
    seed: int = 17,
) -> dict[str, np.ndarray]:
    rng = np.random.default_rng(seed)
    n_modality_dims = n_modality_dims or {
        "transcriptomics": n_kinases, "cnv": n_kinases, "mutation": n_kinases,
        "proteomics": 32, "phosphoproteomics": 192,
    }

    graph = generate_kinase_graph(n_kinases, seed=seed)

    context = rng.integers(0, n_cell_lines, n_samples).astype(np.int64)
    cell_line_programs = rng.normal(0, 0.4, (n_cell_lines, response_dim)).astype(np.float32)

    # each sample perturbs 1 or 2 genes (single or combinatorial), 0 = padding
    combo_sizes = rng.integers(1, max_combo_size + 1, n_samples)
    gene_ids = np.zeros((n_samples, max_combo_size), dtype=np.int64)
    for i, size in enumerate(combo_sizes):
        chosen = rng.choice(np.arange(1, n_genes + 1), size=size, replace=False)  # 1-indexed, 0 reserved for padding
        gene_ids[i, :size] = chosen
    direction_ids = rng.integers(1, n_directions, n_samples).astype(np.int64)
    dose = rng.uniform(0.0, 3.0, n_samples).astype(np.float32)

    # perturbation array used for OOD splitting: primary (first, non-padding) gene id
    perturbation_group = gene_ids[:, 0]

    gene_programs = rng.normal(0, 0.5, (n_genes + 1, response_dim)).astype(np.float32)
    gene_programs[0] = 0.0  # padding contributes nothing
    response = np.zeros((n_samples, response_dim), dtype=np.float32)
    for i in range(n_samples):
        active_genes = gene_ids[i, gene_ids[i] > 0]
        response[i] = gene_programs[active_genes].sum(axis=0) + cell_line_programs[context[i]]
    response += rng.normal(0, 0.15, response.shape).astype(np.float32)

    # redundancy label: proxy for "does the network have slack" -- higher when more
    # OTHER kinases connected to the perturbed gene(s) remain untouched. Replace with
    # real single-vs-combo knockout viability comparisons for biological grounding.
    redundancy_label = np.zeros(n_samples, dtype=np.float32)
    for i in range(n_samples):
        active = set(int(g) for g in gene_ids[i] if g > 0)
        node_idx = [g - 1 for g in active if g - 1 < n_kinases]
        if node_idx:
            neighbor_strength = graph["adjacency"][node_idx].sum(axis=1).mean()
            redundancy_label[i] = 1.0 / (1.0 + np.exp(-2.0 * (neighbor_strength - 1.5)))  # squashed proxy
        else:
            redundancy_label[i] = 0.5

    pathway_projection = rng.normal(size=(response_dim, pathway_dim)).astype(np.float32) * 0.3
    pathway_scores = response @ pathway_projection + rng.normal(0, 0.05, (n_samples, pathway_dim)).astype(np.float32)

    modalities = {
        name: rng.normal(0, 1.0, (n_samples, dim)).astype(np.float32)
        for name, dim in n_modality_dims.items()
    }
    # randomly drop some modalities per sample to exercise the missing-modality masking
    modality_mask = {
        name: (rng.uniform(size=n_samples) > 0.15).astype(np.float32)  # ~15% missing
        for name in n_modality_dims
    }

    data = {
        "context": context,
        "perturbation": perturbation_group,
        "gene_ids": gene_ids,
        "direction_ids": direction_ids,
        "dose": dose,
        "response": response,
        "redundancy_label": redundancy_label,
        "pathway_scores": pathway_scores,
        "kinase_node_features": graph["node_features"],
        "kinase_adjacency": graph["adjacency"],
    }
    for name, arr in modalities.items():
        data[f"modality__{name}"] = arr
    for name, arr in modality_mask.items():
        data[f"modality_mask__{name}"] = arr
    return data


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="data/processed/kinase_synthetic.npz")
    parser.add_argument("--seed", type=int, default=17)
    args = parser.parse_args()
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output, **generate(seed=args.seed))
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
