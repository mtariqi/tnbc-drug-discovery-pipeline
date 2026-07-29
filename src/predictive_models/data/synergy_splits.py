"""
Leakage-safe splitting for drug-PAIR data (synergy), where a single drug identity
can appear in either the drug_row or drug_col position of a given record.

This is NOT a drop-in case for src.predictive_models.data.splits.make_splits(), which assumes
one perturbation identity per row. Naively calling make_splits() on drug_row alone
was checked and found to leak 100% of "held-out" drugs back into training via
drug_col (the same drug, just partnered differently) -- see the project history for
the exact check. This module holds a drug out of BOTH columns simultaneously.

cell_line_id has no such ambiguity (one column, one identity per row), so
unseen_context can still safely reuse src.predictive_models.data.splits.make_splits directly.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class SynergySplit:
    train: np.ndarray
    val: np.ndarray
    test: np.ndarray


def unseen_drug_split(
    drug_row: np.ndarray,
    drug_col: np.ndarray,
    val_fraction: float = 0.15,
    test_fraction: float = 0.20,
    seed: int = 0,
) -> SynergySplit:
    """Partition the UNIQUE SET OF DRUGS (not rows) into train/val/test groups,
    then assign each PAIR-ROW to a split only if BOTH its drugs fall in the same
    group. Rows whose two drugs fall in different groups (e.g. one train-group
    drug paired with one test-group drug) are dropped entirely -- included in
    neither split -- since such a row would either leak a training-set drug's
    identity into a test-set evaluation, or vice versa, if force-assigned to either
    side. This trades away some data for a split that is actually leakage-free,
    rather than approximately so.
    """
    rng = np.random.default_rng(seed)
    all_drugs = np.unique(np.concatenate([drug_row, drug_col]))
    n = len(all_drugs)
    shuffled = rng.permutation(all_drugs)
    n_test = max(1, int(n * test_fraction))
    n_val = max(1, int(n * val_fraction))
    test_drugs = set(shuffled[:n_test].tolist())
    val_drugs = set(shuffled[n_test:n_test + n_val].tolist())
    train_drugs = set(shuffled[n_test + n_val:].tolist())

    def row_group(i: int) -> str | None:
        a, b = int(drug_row[i]), int(drug_col[i])
        if a in train_drugs and b in train_drugs:
            return "train"
        if a in val_drugs and b in val_drugs:
            return "val"
        if a in test_drugs and b in test_drugs:
            return "test"
        return None  # mixed-group pair -- dropped, see docstring

    groups = np.array([row_group(i) for i in range(len(drug_row))], dtype=object)
    train_idx = np.where(groups == "train")[0]
    val_idx = np.where(groups == "val")[0]
    test_idx = np.where(groups == "test")[0]
    dropped = len(drug_row) - len(train_idx) - len(val_idx) - len(test_idx)
    print(f"unseen_drug_split: {len(train_idx)} train / {len(val_idx)} val / {len(test_idx)} test rows; "
          f"{dropped} rows dropped (mixed train/test drug pairs, unassignable without leakage).")

    return SynergySplit(train=train_idx, val=val_idx, test=test_idx)
