from __future__ import annotations

from collections.abc import Mapping

import torch
from torch import nn

from src.predictive_models.encoders import (
    CellStateEncoder,
    KinaseNetworkEncoder,
    PerturbationEncoder,
)


class KinaseCellModel(nn.Module):
    """Fuse cell, perturbation, and RTK/NRTK network representations."""

    def __init__(
        self,
        cell_encoder: CellStateEncoder,
        perturbation_encoder: PerturbationEncoder,
        network_encoder: KinaseNetworkEncoder,
        latent_dim: int,
        response_dim: int,
        pathway_dim: int,
        dropout: float = 0.1,
    ):
        super().__init__()
        self.cell_encoder = cell_encoder
        self.perturbation_encoder = perturbation_encoder
        self.network_encoder = network_encoder
        self.fusion = nn.Sequential(
            nn.LayerNorm(latent_dim * 4),
            nn.Linear(latent_dim * 4, latent_dim * 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(latent_dim * 2, latent_dim),
        )
        self.response_mean = nn.Linear(latent_dim, response_dim)
        self.response_log_variance = nn.Linear(latent_dim, response_dim)
        self.pathway_head = nn.Linear(latent_dim, pathway_dim)
        self.redundancy_head = nn.Sequential(
            nn.Linear(latent_dim, latent_dim // 2),
            nn.GELU(),
            nn.Linear(latent_dim // 2, 1),
        )

    def forward(
        self,
        cell_modalities: Mapping[str, torch.Tensor],
        gene_ids: torch.Tensor,
        direction_ids: torch.Tensor,
        dose: torch.Tensor,
        kinase_node_features: torch.Tensor,
        modality_mask: Mapping[str, torch.Tensor] | None = None,
        gene_features: torch.Tensor | None = None,
        drug_features: torch.Tensor | None = None,
        target_mask: torch.Tensor | None = None,
    ) -> dict[str, torch.Tensor | dict[str, torch.Tensor]]:
        cell = self.cell_encoder(cell_modalities, modality_mask)
        perturbation = self.perturbation_encoder(
            gene_ids, direction_ids, dose, gene_features, drug_features
        )
        network = self.network_encoder(kinase_node_features, target_mask)
        interaction = perturbation.embedding * network.embedding
        fused = self.fusion(
            torch.cat(
                (
                    cell.embedding,
                    perturbation.embedding,
                    network.embedding,
                    interaction,
                ),
                dim=-1,
            )
        )
        return {
            "response_mean": self.response_mean(fused),
            "response_log_variance": self.response_log_variance(fused).clamp(-8, 6),
            "pathway_scores": self.pathway_head(fused),
            "redundancy_score": torch.sigmoid(self.redundancy_head(fused)).squeeze(-1),
            "latent": fused,
            "cell_auxiliary": cell.auxiliary,
            "perturbation_auxiliary": perturbation.auxiliary,
            "network_auxiliary": network.auxiliary,
        }

