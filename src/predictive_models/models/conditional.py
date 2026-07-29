from __future__ import annotations

import math

import torch
from torch import nn


class SinusoidalTimeEmbedding(nn.Module):
    def __init__(self, dim: int):
        super().__init__()
        self.dim = dim

    def forward(self, time: torch.Tensor) -> torch.Tensor:
        half = self.dim // 2
        scale = math.log(10_000) / max(half - 1, 1)
        frequencies = torch.exp(-scale * torch.arange(half, device=time.device))
        embedding = time[:, None] * frequencies[None, :]
        output = torch.cat((embedding.sin(), embedding.cos()), dim=-1)
        return torch.nn.functional.pad(output, (0, self.dim - output.shape[-1]))


class ConditionEncoder(nn.Module):
    def __init__(self, n_contexts: int, n_perturbations: int, c_dim: int, p_dim: int):
        super().__init__()
        self.context = nn.Embedding(n_contexts, c_dim)
        self.perturbation = nn.Embedding(n_perturbations, p_dim)

    def forward(self, context: torch.Tensor, perturbation: torch.Tensor) -> torch.Tensor:
        return torch.cat((self.context(context), self.perturbation(perturbation)), dim=-1)


class ConditionalMLP(nn.Module):
    def __init__(self, n_genes: int, n_contexts: int, n_perturbations: int, hidden_dim: int,
                 context_dim: int, perturbation_dim: int):
        super().__init__()
        self.condition = ConditionEncoder(n_contexts, n_perturbations, context_dim, perturbation_dim)
        total = n_genes + context_dim + perturbation_dim
        self.network = nn.Sequential(
            nn.Linear(total, hidden_dim),
            nn.SiLU(),
            nn.LayerNorm(hidden_dim),
            nn.Linear(hidden_dim, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, n_genes),
        )

    def forward(self, control, context, perturbation):
        condition = self.condition(context, perturbation)
        return self.network(torch.cat((control, condition), dim=-1))


class ConditionalVectorField(nn.Module):
    """Velocity network for conditional flow matching in response-residual space."""

    def __init__(self, n_genes: int, n_contexts: int, n_perturbations: int, hidden_dim: int,
                 context_dim: int, perturbation_dim: int, time_dim: int):
        super().__init__()
        self.condition = ConditionEncoder(n_contexts, n_perturbations, context_dim, perturbation_dim)
        self.time = SinusoidalTimeEmbedding(time_dim)
        total = 2 * n_genes + context_dim + perturbation_dim + time_dim
        self.network = nn.Sequential(
            nn.Linear(total, hidden_dim),
            nn.SiLU(),
            nn.LayerNorm(hidden_dim),
            nn.Linear(hidden_dim, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, n_genes),
        )

    def forward(self, state, time, control, context, perturbation):
        features = torch.cat(
            (state, control, self.condition(context, perturbation), self.time(time)), dim=-1
        )
        return self.network(features)

    def loss(self, response, control, context, perturbation):
        noise = torch.randn_like(response)
        time = torch.rand(response.shape[0], device=response.device)
        state = (1 - time[:, None]) * noise + time[:, None] * response
        target_velocity = response - noise
        return torch.nn.functional.mse_loss(
            self(state, time, control, context, perturbation), target_velocity
        )

    @torch.no_grad()
    def sample(self, control, context, perturbation, steps: int = 100):
        state = torch.randn_like(control)
        dt = 1.0 / steps
        for step in range(steps):
            time = torch.full((state.shape[0],), step / steps, device=state.device)
            state = state + dt * self(state, time, control, context, perturbation)
        return state

