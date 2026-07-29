from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class SplitIndices:
    train: np.ndarray
    val: np.ndarray
    test: np.ndarray


def _partition(values: np.ndarray, val_fraction: float, test_fraction: float, rng):
    unique = np.unique(values).copy()
    rng.shuffle(unique)
    n_test = max(1, round(len(unique) * test_fraction))
    n_val = max(1, round(len(unique) * val_fraction))
    return set(unique[n_test + n_val :]), set(unique[n_test : n_test + n_val]), set(unique[:n_test])


def make_splits(
    context: np.ndarray,
    perturbation: np.ndarray,
    mode: str,
    val_fraction: float,
    test_fraction: float,
    seed: int,
) -> SplitIndices:
    """Create deployment-style splits without group leakage."""
    rng = np.random.default_rng(seed)
    n = len(context)
    if mode == "iid":
        order = rng.permutation(n)
        n_test, n_val = round(n * test_fraction), round(n * val_fraction)
        return SplitIndices(order[n_test + n_val :], order[n_test : n_test + n_val], order[:n_test])

    if mode == "unseen_perturbation":
        train_g, val_g, test_g = _partition(perturbation, val_fraction, test_fraction, rng)
        group = perturbation
        return SplitIndices(
            np.flatnonzero(np.isin(group, list(train_g))),
            np.flatnonzero(np.isin(group, list(val_g))),
            np.flatnonzero(np.isin(group, list(test_g))),
        )

    if mode == "unseen_context":
        train_g, val_g, test_g = _partition(context, val_fraction, test_fraction, rng)
        group = context
        return SplitIndices(
            np.flatnonzero(np.isin(group, list(train_g))),
            np.flatnonzero(np.isin(group, list(val_g))),
            np.flatnonzero(np.isin(group, list(test_g))),
        )

    if mode == "unseen_both":
        c_train, c_val, c_test = _partition(context, val_fraction, test_fraction, rng)
        p_train, p_val, p_test = _partition(perturbation, val_fraction, test_fraction, rng)
        train = np.flatnonzero(np.isin(context, list(c_train)) & np.isin(perturbation, list(p_train)))
        val = np.flatnonzero(np.isin(context, list(c_val)) & np.isin(perturbation, list(p_val)))
        test = np.flatnonzero(np.isin(context, list(c_test)) & np.isin(perturbation, list(p_test)))
        if min(len(train), len(val), len(test)) == 0:
            raise ValueError("Empty split; increase data diversity or adjust fractions")
        return SplitIndices(train, val, test)
    raise ValueError(f"Unknown split mode: {mode}")

