from __future__ import annotations

import torch
from torch import nn

from .common import EncoderOutput, ResidualBlock


def normalize_adjacency(adjacency: torch.Tensor) -> torch.Tensor:
    adjacency = adjacency + torch.eye(
        adjacency.shape[0], device=adjacency.device, dtype=adjacency.dtype
    )
    degree = adjacency.sum(dim=-1).clamp_min(1e-8)
    scale = degree.rsqrt()
    return scale[:, None] * adjacency * scale[None, :]


class GraphConvolution(nn.Module):
    def __init__(self, width: int):
        super().__init__()
        self.projection = nn.Linear(width, width)

    def forward(self, nodes: torch.Tensor, adjacency: torch.Tensor) -> torch.Tensor:
        return torch.einsum("ij,bjd->bid", adjacency, self.projection(nodes))


class KinaseNetworkEncoder(nn.Module):
    """Reason over RTK/NRTK evidence without requiring PyG at inference.

    Node features can contain CTS components, expression, mutation, DepMap,
    survival, druggability, RTK/NRTK class, and ESM-2 summaries. Edge weights can
    combine STRING confidence, co-expression, mutation co-occurrence, PSP substrate
    Jaccard, and the existing bootstrapped redundancy score.
    """

    def __init__(
        self,
        node_feature_dim: int,
        adjacency: torch.Tensor,
        latent_dim: int = 256,
        layers: int = 3,
        dropout: float = 0.1,
    ):
        super().__init__()
        if adjacency.ndim != 2 or adjacency.shape[0] != adjacency.shape[1]:
            raise ValueError("adjacency must be a square [kinase, kinase] matrix")
        self.register_buffer("adjacency", normalize_adjacency(adjacency.float()))
        self.node_projection = nn.Linear(node_feature_dim, latent_dim)
        self.convolutions = nn.ModuleList(
            [GraphConvolution(latent_dim) for _ in range(layers)]
        )
        self.refiners = nn.ModuleList(
            [ResidualBlock(latent_dim, dropout) for _ in range(layers)]
        )
        self.target_attention = nn.Linear(latent_dim, 1)

    def forward(
        self,
        node_features: torch.Tensor,
        target_mask: torch.Tensor | None = None,
    ) -> EncoderOutput:
        if node_features.ndim == 2:
            node_features = node_features.unsqueeze(0)
        if node_features.shape[1] != self.adjacency.shape[0]:
            raise ValueError("node count does not match configured adjacency")
        nodes = self.node_projection(node_features)
        for convolution, refiner in zip(self.convolutions, self.refiners):
            nodes = refiner(nodes + torch.nn.functional.gelu(
                convolution(nodes, self.adjacency)
            ))
        logits = self.target_attention(nodes).squeeze(-1)
        if target_mask is not None:
            logits = logits.masked_fill(
                ~target_mask.bool(), torch.finfo(logits.dtype).min
            )
        weights = torch.softmax(logits, dim=-1)
        graph = (nodes * weights.unsqueeze(-1)).sum(dim=1)
        return EncoderOutput(
            embedding=graph,
            auxiliary={"node_embeddings": nodes, "node_attention": weights},
        )

