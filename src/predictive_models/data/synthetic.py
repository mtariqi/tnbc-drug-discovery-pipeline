from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def generate(
    n_cells: int = 1200,
    n_genes: int = 64,
    n_contexts: int = 6,
    n_perturbations: int = 12,
    seed: int = 17,
) -> dict[str, np.ndarray]:
    """Generate structured perturbation effects with heterogeneous single-cell noise."""
    rng = np.random.default_rng(seed)
    context = rng.integers(0, n_contexts, n_cells)
    perturbation = rng.integers(0, n_perturbations, n_cells)
    context_programs = rng.normal(0, 0.35, (n_contexts, n_genes))
    drug_programs = rng.normal(0, 0.55, (n_perturbations, n_genes))
    interaction = rng.normal(0, 0.18, (n_contexts, n_perturbations, n_genes))

    baseline = 1.5 + context_programs[context] + rng.normal(0, 0.25, (n_cells, n_genes))
    delta = drug_programs[perturbation] + interaction[context, perturbation]
    heteroskedastic = 0.12 + 0.10 * np.abs(delta)
    response = delta + rng.normal(0, heteroskedastic)
    expression = baseline + response

    return {
        "control": baseline.astype(np.float32),
        "response": response.astype(np.float32),
        "expression": expression.astype(np.float32),
        "context": context.astype(np.int64),
        "perturbation": perturbation.astype(np.int64),
        "gene_names": np.asarray([f"GENE_{i:04d}" for i in range(n_genes)]),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="data/processed/synthetic.npz")
    parser.add_argument("--seed", type=int, default=17)
    args = parser.parse_args()
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output, **generate(seed=args.seed))
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()

