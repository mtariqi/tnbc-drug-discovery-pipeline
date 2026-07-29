from __future__ import annotations

from collections.abc import Mapping

import torch
from torch import nn

from .common import EncoderOutput, ResidualBlock


class CellStateEncoder(nn.Module):
    """Encode TNBC multi-omics while explicitly handling missing modalities.

    Inputs are aligned patient, cell-line, or single-cell matrices. Each modality
    is projected independently, then combined with learned modality gates. Missing
    modalities receive a zero mask rather than an imputed biological claim.
    """

    def __init__(
        self,
        modality_dims: Mapping[str, int],
        latent_dim: int = 256,
        dropout: float = 0.1,
    ):
        super().__init__()
        if not modality_dims:
            raise ValueError("At least one cell-state modality is required")
        self.modalities = tuple(sorted(modality_dims))
        self.projectors = nn.ModuleDict(
            {
                name: nn.Sequential(
                    nn.LayerNorm(size),
                    nn.Linear(size, latent_dim),
                    nn.GELU(),
                    nn.Dropout(dropout),
                )
                for name, size in modality_dims.items()
            }
        )
        self.missing_tokens = nn.ParameterDict(
            {name: nn.Parameter(torch.zeros(latent_dim)) for name in self.modalities}
        )
        self.gates = nn.ModuleDict(
            {name: nn.Linear(latent_dim, 1) for name in self.modalities}
        )
        self.refine = ResidualBlock(latent_dim, dropout)

    def forward(
        self,
        modalities: Mapping[str, torch.Tensor],
        modality_mask: Mapping[str, torch.Tensor] | None = None,
    ) -> EncoderOutput:
        if not modalities:
            raise ValueError("No cell-state tensors were supplied")
        batch_size = next(iter(modalities.values())).shape[0]
        encoded, available = [], []
        for name in self.modalities:
            present = name in modalities
            if present:
                state = self.projectors[name](modalities[name])
                mask = (
                    modality_mask[name].float()
                    if modality_mask is not None and name in modality_mask
                    else torch.ones(batch_size, device=state.device)
                )
            else:
                reference = next(iter(modalities.values()))
                state = self.missing_tokens[name].expand(batch_size, -1)
                mask = torch.zeros(batch_size, device=reference.device)
            encoded.append(state)
            available.append(mask)

        states = torch.stack(encoded, dim=1)
        mask = torch.stack(available, dim=1)
        logits = torch.cat(
            [self.gates[name](states[:, index]) for index, name in enumerate(self.modalities)],
            dim=1,
        )
        logits = logits.masked_fill(mask == 0, torch.finfo(logits.dtype).min)
        weights = torch.softmax(logits, dim=1)
        fused = (states * weights.unsqueeze(-1)).sum(dim=1)
        return EncoderOutput(
            embedding=self.refine(fused),
            auxiliary={"modality_weights": weights, "modality_mask": mask},
        )

