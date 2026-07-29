from __future__ import annotations

from dataclasses import dataclass

import torch


@dataclass
class EncoderOutput:
    """Uniform contract shared by every encoder."""

    embedding: torch.Tensor
    auxiliary: dict[str, torch.Tensor]


class ResidualBlock(torch.nn.Module):
    def __init__(self, width: int, dropout: float):
        super().__init__()
        self.block = torch.nn.Sequential(
            torch.nn.LayerNorm(width),
            torch.nn.Linear(width, width * 2),
            torch.nn.GELU(),
            torch.nn.Dropout(dropout),
            torch.nn.Linear(width * 2, width),
        )

    def forward(self, values: torch.Tensor) -> torch.Tensor:
        return values + self.block(values)

