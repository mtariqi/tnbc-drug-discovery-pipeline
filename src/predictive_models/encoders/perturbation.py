from __future__ import annotations

import torch
from torch import nn

from .common import EncoderOutput, ResidualBlock


class PerturbationEncoder(nn.Module):
    """Encode single or combinatorial genetic/drug perturbations.

    Gene IDs are pooled as a permutation-invariant set. Direction distinguishes
    activation, knockout/knockdown, and inhibition; dose is represented on a
    log1p scale. Drug molecular features can be added without changing the API.
    """

    def __init__(
        self,
        n_genes: int,
        n_directions: int = 4,
        gene_feature_dim: int | None = None,
        drug_feature_dim: int | None = None,
        latent_dim: int = 256,
        dropout: float = 0.1,
        padding_idx: int = 0,
    ):
        super().__init__()
        self.padding_idx = padding_idx
        self.gene_embedding = nn.Embedding(
            n_genes, latent_dim, padding_idx=padding_idx
        )
        self.direction_embedding = nn.Embedding(n_directions, latent_dim)
        self.gene_feature_projection = (
            nn.Linear(gene_feature_dim, latent_dim) if gene_feature_dim else None
        )
        self.drug_projection = (
            nn.Linear(drug_feature_dim, latent_dim) if drug_feature_dim else None
        )
        self.dose_projection = nn.Sequential(
            nn.Linear(1, latent_dim), nn.GELU(), nn.Linear(latent_dim, latent_dim)
        )
        self.combine = nn.Sequential(
            nn.Linear(latent_dim * 4, latent_dim),
            nn.GELU(),
            ResidualBlock(latent_dim, dropout),
        )

    def forward(
        self,
        gene_ids: torch.Tensor,
        direction_ids: torch.Tensor,
        dose: torch.Tensor,
        gene_features: torch.Tensor | None = None,
        drug_features: torch.Tensor | None = None,
    ) -> EncoderOutput:
        gene_mask = gene_ids.ne(self.padding_idx)
        gene_states = self.gene_embedding(gene_ids)
        if gene_features is not None:
            if self.gene_feature_projection is None:
                raise ValueError("gene_feature_dim was not configured")
            gene_states = gene_states + self.gene_feature_projection(gene_features)
        denominator = gene_mask.sum(dim=1, keepdim=True).clamp_min(1)
        gene_pool = (gene_states * gene_mask.unsqueeze(-1)).sum(dim=1) / denominator
        direction = self.direction_embedding(direction_ids)
        dose_state = self.dose_projection(torch.log1p(dose.clamp_min(0)).unsqueeze(-1))
        if drug_features is None:
            drug_state = torch.zeros_like(gene_pool)
        else:
            if self.drug_projection is None:
                raise ValueError("drug_feature_dim was not configured")
            drug_state = self.drug_projection(drug_features)
        fused = self.combine(
            torch.cat((gene_pool, direction, dose_state, drug_state), dim=-1)
        )
        return EncoderOutput(
            embedding=fused,
            auxiliary={
                "gene_mask": gene_mask,
                "combination_size": gene_mask.sum(dim=1),
            },
        )

