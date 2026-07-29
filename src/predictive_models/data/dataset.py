from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import torch
from torch.utils.data import Dataset


REQUIRED_KEYS = {"control", "response", "context", "perturbation"}


@dataclass(frozen=True)
class PerturbationArrays:
    control: np.ndarray
    response: np.ndarray
    context: np.ndarray
    perturbation: np.ndarray
    gene_names: np.ndarray

    @classmethod
    def load(cls, path: str | Path) -> "PerturbationArrays":
        data = np.load(path, allow_pickle=False)
        missing = REQUIRED_KEYS.difference(data.files)
        if missing:
            raise ValueError(f"Missing required arrays: {sorted(missing)}")
        n = data["control"].shape[0]
        if data["response"].shape != data["control"].shape:
            raise ValueError("control and response must have identical [cell, gene] shape")
        if data["context"].shape[0] != n or data["perturbation"].shape[0] != n:
            raise ValueError("metadata length does not match cell count")
        genes = data["gene_names"] if "gene_names" in data.files else np.arange(data["control"].shape[1])
        return cls(
            control=data["control"].astype(np.float32),
            response=data["response"].astype(np.float32),
            context=data["context"].astype(np.int64),
            perturbation=data["perturbation"].astype(np.int64),
            gene_names=genes,
        )


class PerturbationDataset(Dataset):
    def __init__(self, arrays: PerturbationArrays, indices: np.ndarray):
        self.arrays = arrays
        self.indices = np.asarray(indices)

    def __len__(self) -> int:
        return len(self.indices)

    def __getitem__(self, item: int) -> tuple[torch.Tensor, ...]:
        idx = self.indices[item]
        return (
            torch.from_numpy(self.arrays.control[idx]),
            torch.tensor(self.arrays.context[idx], dtype=torch.long),
            torch.tensor(self.arrays.perturbation[idx], dtype=torch.long),
            torch.from_numpy(self.arrays.response[idx]),
        )

