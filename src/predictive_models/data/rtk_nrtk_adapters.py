from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


NODE_ALIASES = {
    "gene": ("gene", "gene_symbol", "kinase"),
    "kinase_class": ("kinase_class", "class", "type"),
    "string_centrality": ("string_centrality", "centrality", "ppi_centrality"),
    "depmap_essentiality": ("depmap_essentiality", "essentiality", "crispr_score"),
    "survival_score": ("survival_score", "tcga_survival", "survival"),
    "druggability": ("druggability", "dgidb_score", "drug_score"),
    "cts": ("cts", "composite_target_score", "target_score"),
}

EDGE_ALIASES = {
    "source": ("source", "kinase1", "gene_a", "preferredName_A"),
    "target": ("target", "kinase2", "gene_b", "preferredName_B"),
    "string_score": ("string_score", "combined_score", "score"),
    "coexpression": ("coexpression", "correlation", "pearson_r"),
    "mutation_cooccurrence": ("mutation_cooccurrence", "fisher_score"),
    "substrate_jaccard": ("substrate_jaccard", "jaccard"),
    "redundancy": ("redundancy", "redundancy_score", "composite_redundancy"),
}


def _resolve(frame: pd.DataFrame, aliases: dict[str, tuple[str, ...]]) -> pd.DataFrame:
    lower = {column.lower(): column for column in frame.columns}
    resolved = {}
    for canonical, candidates in aliases.items():
        match = next((lower[item.lower()] for item in candidates if item.lower() in lower), None)
        if match is not None:
            resolved[match] = canonical
    return frame.rename(columns=resolved)


def _robust_scale(values: np.ndarray) -> np.ndarray:
    median = np.nanmedian(values, axis=0, keepdims=True)
    q75 = np.nanpercentile(values, 75, axis=0, keepdims=True)
    q25 = np.nanpercentile(values, 25, axis=0, keepdims=True)
    scaled = (values - median) / np.maximum(q75 - q25, 1e-6)
    return np.nan_to_num(scaled, nan=0.0, posinf=0.0, neginf=0.0)


@dataclass(frozen=True)
class KinaseGraphArrays:
    genes: list[str]
    node_features: np.ndarray
    adjacency: np.ndarray
    node_feature_names: list[str]


def load_kinase_graph(
    node_csv: str | Path,
    edge_csv: str | Path,
    minimum_string_score: float = 700,
) -> KinaseGraphArrays:
    """Adapt current CTS/STRING/PSP/redundancy outputs to encoder tensors."""
    nodes = _resolve(pd.read_csv(node_csv), NODE_ALIASES)
    edges = _resolve(pd.read_csv(edge_csv), EDGE_ALIASES)
    if "gene" not in nodes:
        raise ValueError("Node file needs a gene/gene_symbol/kinase column")
    if not {"source", "target"}.issubset(edges):
        raise ValueError("Edge file needs source/target or kinase1/kinase2 columns")

    genes = nodes["gene"].astype(str).str.upper().tolist()
    index = {gene: position for position, gene in enumerate(genes)}
    numeric_nodes = [
        name for name in (
            "string_centrality", "depmap_essentiality", "survival_score",
            "druggability", "cts"
        ) if name in nodes
    ]
    values = nodes[numeric_nodes].apply(pd.to_numeric, errors="coerce").to_numpy()
    kinase_class = (
        nodes.get("kinase_class", pd.Series("UNKNOWN", index=nodes.index))
        .astype(str).str.upper().eq("RTK").astype(float).to_numpy()[:, None]
    )
    node_features = np.concatenate((_robust_scale(values), kinase_class), axis=1)
    feature_names = numeric_nodes + ["is_rtk"]

    adjacency = np.zeros((len(genes), len(genes)), dtype=np.float32)
    edge_features = [
        name for name in (
            "string_score", "coexpression", "mutation_cooccurrence",
            "substrate_jaccard", "redundancy"
        ) if name in edges
    ]
    for _, edge in edges.iterrows():
        source, target = str(edge["source"]).upper(), str(edge["target"]).upper()
        if source not in index or target not in index:
            continue
        if "string_score" in edges and np.isfinite(pd.to_numeric(edge["string_score"], errors="coerce")):
            if float(edge["string_score"]) < minimum_string_score:
                continue
        evidence = [
            abs(float(pd.to_numeric(edge[name], errors="coerce")))
            for name in edge_features
            if np.isfinite(pd.to_numeric(edge[name], errors="coerce"))
        ]
        weight = float(np.mean(evidence)) if evidence else 1.0
        adjacency[index[source], index[target]] = max(adjacency[index[source], index[target]], weight)
        adjacency[index[target], index[source]] = max(adjacency[index[target], index[source]], weight)
    if adjacency.max() > 0:
        adjacency /= adjacency.max()
    return KinaseGraphArrays(genes, node_features.astype(np.float32), adjacency, feature_names)

