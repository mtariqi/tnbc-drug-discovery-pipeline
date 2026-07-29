"""
SynergyGNN: predicts real drug-pair synergy metrics (Bliss, Loewe, HSA, ZIP, per
DrugComb's convention) from two drugs' molecular graphs and the cell-line context
they were tested in.

This is the "true predictive model" replacing PairCTS/TripletCTS's current
heuristic synergy placeholder (explicitly named in the manuscript as "the most
consequential planned extension"). It does NOT replace CTS/PairCTS/TripletCTS
themselves -- the intended integration is: once trained and validated, this
model's predicted synergy score becomes the real "synergy evidence" term HCOS
currently has no data for, rather than a wholesale replacement of the composite
scoring approach.

KNOWN LIMITATION, stated plainly rather than glossed over: this model only handles
drugs with a valid small-molecule SMILES string. Biologics -- including
trastuzumab, which appears in this project's own headline focal-patient regimen --
have no meaningful SMILES and cannot be featurized this way. Any regimen
containing a biologic needs either (a) exclusion from this model's real-data
training/evaluation, or (b) a separate, non-graph representation for biologics
(e.g. a learned per-drug ID embedding, matching how PerturbationEncoder already
handles genes) merged in as an alternative input path. This is not yet implemented
-- decide before training on your real DrugComb pull whether trastuzumab-containing
rows should be dropped or handled via a fallback embedding.
"""
from __future__ import annotations

import torch
import torch.nn as nn

from src.predictive_models.encoders.molecule import MoleculeGNNEncoder


class SynergyGNN(nn.Module):
    def __init__(
        self,
        atom_feature_dim: int,
        n_cell_lines: int,
        latent_dim: int = 128,
        gnn_layers: int = 3,
        cell_embedding_dim: int = 32,
        synergy_metrics: tuple[str, ...] = ("bliss", "loewe", "hsa", "zip"),
        dropout: float = 0.1,
    ):
        super().__init__()
        # Shared weights across drug_row / drug_col: synergy is a property of the
        # PAIR, and there's no principled reason drug identity should be encoded
        # differently depending on which column it happened to land in.
        self.molecule_encoder = MoleculeGNNEncoder(atom_feature_dim, latent_dim, gnn_layers, dropout)

        # v1 cell-line representation: a learned ID embedding, same pattern as
        # PerturbationEncoder's gene embedding. Swap for real DepMap expression/CNV
        # features (reusing CellStateEncoder) once cell-line-level omics are wired
        # in here -- this is a deliberately simple starting point, not a final design.
        self.cell_embedding = nn.Embedding(n_cell_lines, cell_embedding_dim)

        self.synergy_metrics = synergy_metrics
        # Symmetric combination -- synergy(A,B) must equal synergy(B,A), since it's a
        # property of the unordered pair, not drug_row vs. drug_col specifically
        # (that ordering is an artifact of how DrugComb happens to store each row).
        # A naive concat([row, col, ...]) would give DIFFERENT predictions for
        # (A,B) vs (B,A) -- caught and fixed here before this went any further.
        # sum(row,col) and |row-col| are both genuinely order-invariant; interaction
        # (elementwise product) is too.
        combined_dim = latent_dim * 3 + cell_embedding_dim  # sum + abs_diff + interaction + cell
        self.head = nn.Sequential(
            nn.Linear(combined_dim, latent_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(latent_dim, len(synergy_metrics)),
        )

    def forward(
        self,
        drug_row_atoms: torch.Tensor, drug_row_adjacency: torch.Tensor, drug_row_mask: torch.Tensor,
        drug_col_atoms: torch.Tensor, drug_col_adjacency: torch.Tensor, drug_col_mask: torch.Tensor,
        cell_line_id: torch.Tensor,
    ) -> torch.Tensor:
        row_embedding = self.molecule_encoder(drug_row_atoms, drug_row_adjacency, drug_row_mask)
        col_embedding = self.molecule_encoder(drug_col_atoms, drug_col_adjacency, drug_col_mask)
        pair_sum = row_embedding + col_embedding
        pair_absdiff = (row_embedding - col_embedding).abs()
        interaction = row_embedding * col_embedding
        cell_embedding = self.cell_embedding(cell_line_id)

        combined = torch.cat([pair_sum, pair_absdiff, interaction, cell_embedding], dim=-1)
        return self.head(combined)  # [B, len(synergy_metrics)]
