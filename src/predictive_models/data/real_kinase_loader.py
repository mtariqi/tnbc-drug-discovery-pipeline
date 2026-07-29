"""
Real-data loader: wires depmap_multiomic_loader.py / depmap_supplemental_loader.py
(your own code, confirmed against real DepMap files this session) into the exact
array contract src/predictive_models/train_kinase.py expects -- the same keys
data/kinase_synthetic.py produces, so train_kinase.py itself needs no changes.

Drop this file alongside depmap_multiomic_loader.py / depmap_supplemental_loader.py
(e.g. in your real project's src/data_loaders/, or anywhere on PYTHONPATH) since it
imports them directly, unmodified.

WHAT THIS BUILDS FROM REAL DATA (confirmed against your uploaded loader code):
  - response                   <- CRISPRGeneDependency.csv, parse_depmap_dependency()
                                   (a real scalar per (cell_line, kinase) row --
                                   matches configs/kinase_cell.yaml's response_dim=1)
  - modality__transcriptomics   <- OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv
  - modality__cnv                <- OmicsCNGeneWGS.csv
  - modality__mutation            <- OmicsSomaticMutationsMatrixDamaging.csv
  - context = cell_line_id (ModelID, int-encoded), perturbation = kinase_id
    (gene symbol, int-encoded, 1-indexed so 0 stays the padding index)
  all restricted to your real 90-kinase panel and confirmed-TNBC cell-line list.

WHAT THIS DOES NOT BUILD, AND WHY -- STATED PLAINLY, NOT GUESSED AT:
  - kinase_node_features / kinase_adjacency: these belong to your already-validated
    CTS pipeline (kinase_scoring_pipeline.py / build_and_verify_real_cts.py) and the
    cached STRING edges file. I have not been given those two files' exact real
    column layouts this session -- only their existence and rough description
    (`CTS_all_90_kinases.tsv`, `string_edges.tsv`, 464 edges). Rather than guess
    column names, `load_kinase_graph_from_cts()` below follows the same discipline
    your own `inspect_depmap_csv()` uses: call it with `inspect_only=True` first and
    send me what it prints if the column names don't match its first guess -- don't
    let it silently load the wrong column.
  - redundancy_label / pathway_scores: correctly absent. No real label exists yet
    per docs/limitations.md, and configs/kinase_cell.yaml already sets both
    loss_weights to 0, so train_kinase.py doesn't require them.
  - proteomics / phosphoproteomics modalities: correctly all-missing for every row
    produced here. DepMap cell lines and CPTAC patients are different populations
    with no shared identifier -- passed as zero-filled arrays with an all-zero
    modality_mask, which CellStateEncoder's missing-modality gating already handles
    (see src/predictive_models/encoders/cell_state.py) rather than needing new code for this.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from depmap_multiomic_loader import (
    inspect_depmap_csv,
    parse_depmap_cnv,
    parse_depmap_damaging_mutations,
    parse_depmap_dependency,
    parse_depmap_expression,
)
from depmap_supplemental_loader import flag_pan_essential_kinases


def load_kinase_panel(path: str | Path) -> list[str]:
    with open(path) as handle:
        return [line.strip() for line in handle if line.strip()]


def load_tnbc_cell_lines(path: str | Path) -> list[str]:
    with open(path) as handle:
        return [line.strip() for line in handle if line.strip()]


def build_real_response_and_modalities(
    kinase_panel: list[str],
    tnbc_cell_lines: list[str],
    expression_csv: str,
    cnv_csv: str,
    dependency_csv: str,
    damaging_mutation_csv: str,
    knockout_direction_id: int,
) -> dict[str, np.ndarray]:
    """One row per (cell_line, kinase) pair -- exactly dependency_ml_model.py's
    training contract, reused here rather than re-derived, since it's already
    correct: grouped by cell_line_id, dependency_prob as the real label.
    """
    print("Loading real DepMap files -- this can take a minute on the full "
          "pan-cancer CSVs before TNBC/kinase-panel restriction.")
    expr = parse_depmap_expression(expression_csv, kinase_panel, cell_line_ids=tnbc_cell_lines)
    cnv = parse_depmap_cnv(cnv_csv, kinase_panel, cell_line_ids=tnbc_cell_lines)
    dep = parse_depmap_dependency(dependency_csv, kinase_panel, cell_line_ids=tnbc_cell_lines)
    mut = parse_depmap_damaging_mutations(damaging_mutation_csv, kinase_panel, cell_line_ids=tnbc_cell_lines)

    # Inner-join on (cell_line_id, kinase_id): a row only exists if the label
    # (dependency_prob) AND all three cell-state features are present -- matches
    # dependency_ml_model.py's "drop rows with no label" policy for the target,
    # but note this is stricter than that function on the FEATURE side (it
    # left-joins + median-imputes missing features; doing a real inner join here
    # instead so every real value reported is an actually-measured value, not an
    # imputed placeholder -- revisit if this drops too many rows in practice).
    merged = (
        dep[["cell_line_id", "kinase_id", "dependency_prob"]]
        .merge(expr[["cell_line_id", "kinase_id", "log_tpm"]], on=["cell_line_id", "kinase_id"])
        .merge(cnv[["cell_line_id", "kinase_id", "log2_cn"]], on=["cell_line_id", "kinase_id"])
        .merge(mut[["cell_line_id", "kinase_id", "damaging_mutation"]], on=["cell_line_id", "kinase_id"])
    )
    print(f"{len(merged)} real (cell_line, kinase) rows with a complete "
          f"dependency_prob + expression + CNV + mutation record "
          f"(out of {len(tnbc_cell_lines)} cell lines x {len(kinase_panel)} kinases = "
          f"{len(tnbc_cell_lines) * len(kinase_panel)} possible).")

    cell_lines_present = sorted(merged["cell_line_id"].unique())
    cell_line_to_idx = {c: i for i, c in enumerate(cell_lines_present)}
    kinase_to_idx = {k: i + 1 for i, k in enumerate(kinase_panel)}  # 1-indexed, 0 = padding

    context = merged["cell_line_id"].map(cell_line_to_idx).to_numpy(dtype=np.int64)
    perturbation = merged["kinase_id"].map(kinase_to_idx).to_numpy(dtype=np.int64)
    gene_ids = np.stack([perturbation, np.zeros_like(perturbation)], axis=1)  # single knockout, no combos
    direction_ids = np.full(len(merged), knockout_direction_id, dtype=np.int64)
    dose = np.zeros(len(merged), dtype=np.float32)  # CRISPR knockout has no pharmacological dose
    response = merged["dependency_prob"].to_numpy(dtype=np.float32)[:, None]  # shape (n, 1)

    # Per-cell-line profile ACROSS THE WHOLE PANEL (not just the perturbed gene) --
    # this is CellStateEncoder's "cell state" input, distinct from the specific
    # gene being perturbed, which the perturbation/network encoders handle.
    def pivot_to_panel_profile(df: pd.DataFrame, value_col: str) -> tuple[np.ndarray, np.ndarray]:
        wide = df.pivot(index="cell_line_id", columns="kinase_id", values=value_col)
        wide = wide.reindex(index=cell_lines_present, columns=kinase_panel)
        n_missing = int(wide.isna().sum().sum())
        if n_missing:
            column_median = wide.median()
            # Same class of bug already found and fixed in this project's own Track C
            # work: an entirely-missing column's median() is NaN, not a fillable
            # value -- naive fillna(median()) would silently leave it as NaN. Catch
            # that case explicitly rather than reproducing it here.
            entirely_missing = column_median[column_median.isna()].index.tolist()
            if entirely_missing:
                print(f"WARNING: {len(entirely_missing)} kinase(s) have zero '{value_col}' "
                      f"coverage across every cell line in this cohort: {entirely_missing}. "
                      f"Filling with 0.0 as an explicit no-data placeholder (not a real "
                      f"imputed value) -- treat these kinases' {value_col} feature as "
                      f"uninformative for this modality, not as a real zero measurement.")
                column_median = column_median.fillna(0.0)
            wide = wide.fillna(column_median)
            print(f"Imputed {n_missing} missing '{value_col}' panel values with the column median "
                  f"(matches dependency_ml_model.py's existing feature-imputation policy).")
        profile = wide.loc[merged["cell_line_id"]].to_numpy(dtype=np.float32)
        mask = np.ones(len(merged), dtype=np.float32)
        return profile, mask

    transcriptomics, transcriptomics_mask = pivot_to_panel_profile(expr, "log_tpm")
    cnv_profile, cnv_mask = pivot_to_panel_profile(cnv, "log2_cn")
    mutation_profile, mutation_mask = pivot_to_panel_profile(mut, "damaging_mutation")

    n = len(merged)
    data = {
        "context": context,
        "perturbation": perturbation,
        "gene_ids": gene_ids,
        "direction_ids": direction_ids,
        "dose": dose,
        "response": response,
        "modality__transcriptomics": transcriptomics,
        "modality_mask__transcriptomics": transcriptomics_mask,
        "modality__cnv": cnv_profile,
        "modality_mask__cnv": cnv_mask,
        "modality__mutation": mutation_profile,
        "modality_mask__mutation": mutation_mask,
        # Correctly all-missing -- see module docstring. Shapes match
        # configs/kinase_cell.yaml's proteomics=32 / phosphoproteomics=192.
        "modality__proteomics": np.zeros((n, 32), dtype=np.float32),
        "modality_mask__proteomics": np.zeros(n, dtype=np.float32),
        "modality__phosphoproteomics": np.zeros((n, 192), dtype=np.float32),
        "modality_mask__phosphoproteomics": np.zeros(n, dtype=np.float32),
    }
    return data


def load_kinase_graph_from_cts(
    cts_tsv_path: str | Path,
    string_edges_tsv_path: str | Path,
    kinase_panel: list[str],
    inspect_only: bool = False,
) -> dict[str, np.ndarray] | None:
    """UNVERIFIED COLUMN NAMES -- inspect before trusting.

    Best-guess column names based on the project description
    ("CTS_all_90_kinases.tsv", "string_edges.tsv", 464 edges, centrality/
    essentiality/survival/druggability/cts columns) -- NOT confirmed against the
    real files the way depmap_multiomic_loader.py's columns were. Run with
    inspect_only=True first; if the printed columns don't match the guesses in
    this function's body, tell me what they actually are rather than letting this
    silently pick the wrong column, the same mistake your own project's
    docs/limitations.md flagged and fixed for the mutation-matrix format.
    """
    cts_df = pd.read_csv(cts_tsv_path, sep="\t")
    edges_df = pd.read_csv(string_edges_tsv_path, sep="\t")

    if inspect_only:
        print(f"=== {cts_tsv_path} ===\ncolumns: {cts_df.columns.tolist()}\n{cts_df.head(3)}\n")
        print(f"=== {string_edges_tsv_path} ===\ncolumns: {edges_df.columns.tolist()}\n{edges_df.head(3)}\n")
        return None

    # Best guesses -- CONFIRM against the printed output above before trusting this.
    node_cols = ["string_centrality", "depmap_essentiality", "survival_score", "dgidb_score", "cts", "is_rtk"]
    missing_cols = [c for c in node_cols if c not in cts_df.columns]
    if missing_cols:
        raise ValueError(
            f"Guessed columns {missing_cols} not found in {cts_tsv_path}. "
            f"Real columns are: {cts_df.columns.tolist()}. Run with inspect_only=True "
            f"and adjust node_cols above to match -- do not guess further."
        )

    cts_df = cts_df.set_index(cts_df.columns[0]).reindex(kinase_panel)
    node_features = cts_df[node_cols].to_numpy(dtype=np.float32)

    n = len(kinase_panel)
    idx = {g: i for i, g in enumerate(kinase_panel)}
    adjacency = np.zeros((n, n), dtype=np.float32)
    source_col, target_col, weight_col = edges_df.columns[:3]  # CONFIRM against inspect_only output
    for _, row in edges_df.iterrows():
        s, t = row[source_col], row[target_col]
        if s in idx and t in idx:
            w = float(row[weight_col])
            adjacency[idx[s], idx[t]] = max(adjacency[idx[s], idx[t]], w)
            adjacency[idx[t], idx[s]] = max(adjacency[idx[t], idx[s]], w)
    if adjacency.max() > 0:
        adjacency /= adjacency.max()

    return {"kinase_node_features": node_features, "kinase_adjacency": adjacency}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", default=str(Path.home() / "rtk_nrtk_tnbc"),
                         help="Real data root, matching run_depmap.py's BASE convention.")
    parser.add_argument("--output", default="data/processed/kinase_real.npz")
    parser.add_argument("--inspect-cts-only", action="store_true",
                         help="Print CTS/STRING file columns and exit, without building arrays.")
    args = parser.parse_args()
    base = Path(args.base)

    kinase_panel = load_kinase_panel(base / "data/raw/kinases/kinase_90_list.txt")
    tnbc_cell_lines = load_tnbc_cell_lines(base / "data/raw/depmap/tnbc_model_ids.txt")
    print(f"{len(kinase_panel)}-kinase panel, {len(tnbc_cell_lines)} confirmed TNBC cell lines.")

    if args.inspect_cts_only:
        load_kinase_graph_from_cts(
            base / "data/processed/scoring/CTS_all_90_kinases.tsv",  # CONFIRM real path
            base / "data/processed/string/string_edges.tsv",
            kinase_panel, inspect_only=True,
        )
        return

    data = build_real_response_and_modalities(
        kinase_panel=kinase_panel,
        tnbc_cell_lines=tnbc_cell_lines,
        expression_csv=str(base / "data/raw/depmap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"),
        cnv_csv=str(base / "data/raw/depmap/OmicsCNGeneWGS.csv"),
        dependency_csv=str(base / "data/raw/depmap/CRISPRGeneDependency.csv"),
        damaging_mutation_csv=str(base / "data/raw/depmap/OmicsSomaticMutationsMatrixDamaging.csv"),
        knockout_direction_id=2,  # matches configs/kinase_cell.yaml perturbation.directions.knockout_or_knockdown
    )

    graph = load_kinase_graph_from_cts(
        base / "data/processed/scoring/CTS_all_90_kinases.tsv",
        base / "data/processed/string/string_edges.tsv",
        kinase_panel,
    )
    data.update(graph)

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output, **data)
    print(f"Wrote {output} -- point configs/kinase_cell.yaml's data.path at this file.")


if __name__ == "__main__":
    main()
