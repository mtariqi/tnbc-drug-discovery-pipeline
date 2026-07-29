"""
SMILES -> fixed-size molecular graph featurization.

Runs offline (RDKit is CPU-only and not torch-native), producing cached
atom-feature / adjacency / atom-mask arrays per unique drug, stored once in the
processed dataset rather than recomputed per training step. This mirrors how
kinase_node_features/kinase_adjacency are precomputed once for the whole panel in
the KinaseCell pipeline, rather than rebuilt per batch.

Design choice: fixed max_atoms padding (default 60) rather than variable-size graphs
with a batching library (e.g. torch_geometric's Batch). This keeps the molecule GNN
in the same plain-PyTorch style as KinaseNetworkEncoder (explicitly built to avoid a
PyG dependency at inference -- see encoders/kinase_network.py's docstring), at the
cost of wasted compute on padding atoms for small molecules. 60 atoms comfortably
covers typical small-molecule drugs (afatinib, alpelisib, trastuzumab is a biologic
and won't featurize this way -- see note below); revisit if your real DrugComb
compound set includes larger macrocycles or biologics.

KNOWN LIMITATION: this featurizer only handles small-molecule drugs with a valid
SMILES string. Biologics (e.g. trastuzumab, a monoclonal antibody) have no
meaningful small-molecule SMILES and will fail here -- they need a different
representation entirely (e.g. a learned ID embedding, same as the gene/kinase
embeddings elsewhere in this codebase) or must be excluded from the synergy model's
training set. Check which of your real DrugComb drugs are biologics before running
this at scale.
"""
from __future__ import annotations

import numpy as np
import torch
import torch.nn as nn

from src.predictive_models.encoders.common import ResidualBlock

try:
    from rdkit import Chem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

# Atom feature vocabulary -- deliberately small and interpretable rather than an
# exhaustive descriptor set, matching the level of complexity already used
# elsewhere in this codebase (e.g. the 6-feature kinase node vector).
ATOM_VOCAB = ["C", "N", "O", "S", "F", "Cl", "Br", "I", "P", "B"]  # common drug-like elements
HYBRIDIZATION_VOCAB = ["SP", "SP2", "SP3", "SP3D", "SP3D2"]


def atom_features(atom) -> np.ndarray:
    """One atom -> a fixed-length feature vector: one-hot element (+ 'other'),
    degree, formal charge, aromaticity, one-hot hybridization (+ 'other'), in-ring.
    """
    symbol = atom.GetSymbol()
    element_onehot = [1.0 if symbol == s else 0.0 for s in ATOM_VOCAB]
    element_onehot.append(1.0 if symbol not in ATOM_VOCAB else 0.0)  # 'other' bucket

    hyb = str(atom.GetHybridization())
    hyb_onehot = [1.0 if hyb == h else 0.0 for h in HYBRIDIZATION_VOCAB]
    hyb_onehot.append(1.0 if hyb not in HYBRIDIZATION_VOCAB else 0.0)

    scalar_features = [
        float(atom.GetDegree()) / 4.0,          # normalize roughly to [0,1]
        float(atom.GetFormalCharge()),
        1.0 if atom.GetIsAromatic() else 0.0,
        1.0 if atom.IsInRing() else 0.0,
    ]
    return np.array(element_onehot + hyb_onehot + scalar_features, dtype=np.float32)


ATOM_FEATURE_DIM = len(ATOM_VOCAB) + 1 + len(HYBRIDIZATION_VOCAB) + 1 + 4  # = 21


def featurize_smiles(smiles: str, max_atoms: int = 60) -> dict[str, np.ndarray] | None:
    """Returns {"atom_features": [max_atoms, ATOM_FEATURE_DIM], "adjacency":
    [max_atoms, max_atoms], "atom_mask": [max_atoms]} or None if the SMILES is
    invalid or the molecule exceeds max_atoms (caller should log and exclude
    such drugs rather than silently truncate the graph).
    """
    if not RDKIT_AVAILABLE:
        raise ImportError(
            "RDKit is required for SMILES featurization (pip install rdkit). "
            "Not available in this environment -- this function was written but "
            "could not be executed here to verify against a real molecule."
        )
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    n_atoms = mol.GetNumAtoms()
    if n_atoms > max_atoms:
        return None

    feats = np.zeros((max_atoms, ATOM_FEATURE_DIM), dtype=np.float32)
    adjacency = np.zeros((max_atoms, max_atoms), dtype=np.float32)
    mask = np.zeros(max_atoms, dtype=np.float32)

    for i, atom in enumerate(mol.GetAtoms()):
        feats[i] = atom_features(atom)
        mask[i] = 1.0

    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        # bond order as edge weight (single=1, double=2, triple=3, aromatic=1.5) --
        # a simple, interpretable choice; swap for a learned bond-type embedding if
        # this turns out to matter empirically.
        order = bond.GetBondTypeAsDouble()
        adjacency[i, j] = order
        adjacency[j, i] = order

    return {"atom_features": feats, "adjacency": adjacency, "atom_mask": mask}


def featurize_drug_table(smiles_by_drug: dict[str, str], max_atoms: int = 60) -> dict[str, dict]:
    """Featurize a {drug_name: smiles} table, returning only drugs that
    featurized successfully. Prints which drugs failed and why (invalid SMILES vs.
    too large) so failures are visible rather than silently dropped -- the same
    discipline as the missing-gene/missing-modality warnings elsewhere in this
    codebase.
    """
    results = {}
    failed_invalid, failed_too_large = [], []
    for name, smiles in smiles_by_drug.items():
        if not RDKIT_AVAILABLE:
            raise ImportError("RDKit required -- see featurize_smiles docstring.")
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            failed_invalid.append(name)
            continue
        if mol.GetNumAtoms() > max_atoms:
            failed_too_large.append(name)
            continue
        results[name] = featurize_smiles(smiles, max_atoms=max_atoms)

    if failed_invalid:
        print(f"Warning: {len(failed_invalid)} drug(s) had no valid SMILES (likely "
              f"biologics, e.g. antibodies like trastuzumab -- these need a different "
              f"representation, not this featurizer): {failed_invalid}")
    if failed_too_large:
        print(f"Warning: {len(failed_too_large)} drug(s) exceeded max_atoms={max_atoms}, "
              f"excluded rather than truncated: {failed_too_large}")
    return results


class MoleculeGNNEncoder(nn.Module):
    """Mirrors KinaseNetworkEncoder's design (symmetric-normalized adjacency
    propagation + masked attention pooling) for architectural consistency, and for
    the same reason that encoder avoids torch_geometric: fewer moving dependency
    parts at inference time.

    Unlike KinaseNetworkEncoder (one FIXED graph repeated across a batch), each
    molecule here has its OWN graph, already naturally batched -- adjacency shape
    is [batch, max_atoms, max_atoms] from the start, no repeat/expand trick needed.
    """

    def __init__(self, atom_feature_dim: int, latent_dim: int, layers: int = 3, dropout: float = 0.1):
        super().__init__()
        self.node_projection = nn.Linear(atom_feature_dim, latent_dim)
        self.convolutions = nn.ModuleList([nn.Linear(latent_dim, latent_dim) for _ in range(layers)])
        self.refiners = nn.ModuleList([ResidualBlock(latent_dim, dropout) for _ in range(layers)])
        self.attention = nn.Linear(latent_dim, 1)

    def _normalize_adjacency(self, adjacency: torch.Tensor) -> torch.Tensor:
        """Symmetric normalization D^-1/2 A D^-1/2, per molecule in the batch, with
        self-loops added so an atom's own features survive propagation."""
        batch_size, n, _ = adjacency.shape
        eye = torch.eye(n, device=adjacency.device).unsqueeze(0).expand(batch_size, -1, -1)
        a_hat = adjacency + eye
        degree = a_hat.sum(dim=-1)
        d_inv_sqrt = torch.where(degree > 0, degree.pow(-0.5), torch.zeros_like(degree))
        d_mat = torch.diag_embed(d_inv_sqrt)
        return torch.bmm(torch.bmm(d_mat, a_hat), d_mat)

    def forward(self, atom_features: torch.Tensor, adjacency: torch.Tensor, atom_mask: torch.Tensor) -> torch.Tensor:
        """atom_features: [B, max_atoms, F], adjacency: [B, max_atoms, max_atoms],
        atom_mask: [B, max_atoms] (1 = real atom, 0 = padding). Returns [B, latent_dim].
        """
        norm_adj = self._normalize_adjacency(adjacency)
        nodes = self.node_projection(atom_features)
        for conv, refine in zip(self.convolutions, self.refiners):
            propagated = torch.bmm(norm_adj, conv(nodes))
            nodes = refine(torch.relu(propagated))

        logits = self.attention(nodes).squeeze(-1)  # [B, max_atoms]
        logits = logits.masked_fill(atom_mask == 0, torch.finfo(logits.dtype).min)
        weights = torch.softmax(logits, dim=-1)
        molecule_embedding = (nodes * weights.unsqueeze(-1)).sum(dim=1)
        return molecule_embedding
