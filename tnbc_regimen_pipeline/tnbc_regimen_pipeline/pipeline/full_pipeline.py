"""
Full pipeline orchestration.

Wires tnbc_regimen_pipeline.discovery and tnbc_regimen_pipeline.cohort
together: parallel literature discovery for coverage-gap genes ->
provenance-tracked merge -> cohort-wide MDCOE/HCOS scoring -> frequency
aggregation -> (optional) before/after diff -> (optional) export ->
(optional) visual summary. Neither the discovery nor cohort packages are
modified by anything here; this module only composes them.
"""

from __future__ import annotations

import json
import logging
import os
from datetime import datetime, timezone
from typing import Callable, Dict, List, Optional

import pandas as pd

from ..discovery.parallelization import parallel_discovery
from ..discovery.provenance import build_merged_pool_with_provenance
from ..cohort.cohort_wide_regimen_analysis import run_cohort_wide_analysis, aggregate_regimen_frequency
from ..utils.validation import validate_gene_drug_schema
from ..config import PipelineConfig
from .diff import diff_regimens_with_and_without_discovery
from .reproducibility import PIPELINE_VERSION, get_git_commit, stamp_dataframe

logger = logging.getLogger(__name__)


# =====================================================================
# EXPORTS
# =====================================================================

def export_results(
    output_dir: str,
    discovery_result: pd.DataFrame,
    merged_gene_drugs: Dict[str, List[str]],
    provenance: Dict[str, Dict[str, List[str]]],
    cohort_result: pd.DataFrame,
    freq: pd.DataFrame,
    diff: Optional[pd.DataFrame] = None,
) -> Dict[str, str]:
    """Writes every intermediate/final artifact as CSV (tabular) or JSON
    (nested), all stamped with pipeline version/git commit/timestamp.
    Returns a dict of {name: path} for what was written."""
    os.makedirs(output_dir, exist_ok=True)
    paths = {}

    def _write_csv(df: pd.DataFrame, name: str):
        p = os.path.join(output_dir, f"{name}.csv")
        stamp_dataframe(df).to_csv(p, index=False)
        paths[name] = p

    _write_csv(discovery_result, "discovery_result")
    _write_csv(cohort_result, "cohort_result")
    _write_csv(freq, "regimen_frequency")
    if diff is not None:
        _write_csv(diff, "regimen_diff")

    meta = {
        "pipeline_version": PIPELINE_VERSION,
        "git_commit": get_git_commit(),
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
    }
    for name, obj in [("merged_gene_drugs", merged_gene_drugs), ("provenance", provenance)]:
        p = os.path.join(output_dir, f"{name}.json")
        with open(p, "w") as f:
            json.dump({"_metadata": meta, "data": obj}, f, indent=2)
        paths[name] = p

    logger.info(f"exported {len(paths)} artifact(s) to {output_dir}")
    return paths


# =====================================================================
# VISUAL SUMMARY
# =====================================================================

def generate_frequency_bar_chart(freq: pd.DataFrame, output_path: str) -> str:
    """Horizontal bar chart of regimen -> #patients top-ranked. Uses
    matplotlib's non-interactive 'Agg' backend (no display in typical
    pipeline-run environments)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, max(3, 0.45 * len(freq))))
    if freq.empty:
        ax.text(0.5, 0.5, "No regimens to display", ha="center", va="center")
    else:
        ordered = freq.sort_values("n_patients_top_ranked")
        ax.barh(ordered["regimen"], ordered["n_patients_top_ranked"])
        ax.set_xlabel("Patients with this top-ranked regimen")
        ax.set_title("Regimen frequency across cohort")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    return output_path


def generate_gene_drug_regimen_flow_svg(
    provenance: Dict[str, Dict[str, List[str]]],
    cohort_result: pd.DataFrame,
    output_path: str,
) -> str:
    """
    Hand-rolled 3-column flow diagram (gene -> drug -> regimen). This is a
    simplification, not a true Sankey (plotly was not available when this
    was built): gene->drug edges just mark "contributes to"; drug->regimen
    edges are weighted by how many patients' top regimen contained that
    drug. Good for seeing which genes fed which drugs into which
    regimens; not a substitute for a real Sankey if that level of rigor
    matters.
    """
    genes = sorted(provenance.keys())
    drugs = sorted({drug for gene_drugs in provenance.values() for drug in gene_drugs})

    valid_regimens = cohort_result.dropna(subset=["top_regimen"])
    drug_to_regimen_count: Dict[str, int] = {}
    for regimen in valid_regimens["top_regimen"]:
        for drug in regimen.split(" + "):
            drug_to_regimen_count[drug] = drug_to_regimen_count.get(drug, 0) + 1
    regimens = sorted(valid_regimens["top_regimen"].unique())

    width, height = 900, max(300, 40 * max(len(genes), len(drugs), len(regimens), 1))
    col_x = {"gene": 80, "drug": width // 2, "regimen": width - 80}

    def _y_positions(items):
        n = len(items)
        if n == 0:
            return {}
        step = height / (n + 1)
        return {item: step * (i + 1) for i, item in enumerate(items)}

    gene_y = _y_positions(genes)
    drug_y = _y_positions(drugs)
    regimen_y = _y_positions(regimens)

    svg_parts = [f'<svg viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg" font-family="sans-serif" font-size="12">']
    svg_parts.append(f'<rect width="{width}" height="{height}" fill="white"/>')

    for gene, gene_drugs in provenance.items():
        if gene not in gene_y:
            continue
        for drug in gene_drugs:
            if drug not in drug_y:
                continue
            literature_sourced = len(gene_drugs[drug]) > 0
            stroke = "#c0392b" if literature_sourced else "#95a5a6"
            svg_parts.append(
                f'<line x1="{col_x["gene"]}" y1="{gene_y[gene]:.1f}" '
                f'x2="{col_x["drug"]}" y2="{drug_y[drug]:.1f}" '
                f'stroke="{stroke}" stroke-width="1.5" opacity="0.6"/>'
            )

    for regimen in regimens:
        for drug in regimen.split(" + "):
            if drug not in drug_y or regimen not in regimen_y:
                continue
            weight = drug_to_regimen_count.get(drug, 1)
            svg_parts.append(
                f'<line x1="{col_x["drug"]}" y1="{drug_y[drug]:.1f}" '
                f'x2="{col_x["regimen"]}" y2="{regimen_y[regimen]:.1f}" '
                f'stroke="#2980b9" stroke-width="{1 + weight}" opacity="0.5"/>'
            )

    for label, y in gene_y.items():
        svg_parts.append(f'<circle cx="{col_x["gene"]}" cy="{y:.1f}" r="4" fill="#2c3e50"/>')
        svg_parts.append(f'<text x="{col_x["gene"] - 10}" y="{y:.1f}" text-anchor="end" dominant-baseline="middle">{label}</text>')
    for label, y in drug_y.items():
        svg_parts.append(f'<circle cx="{col_x["drug"]}" cy="{y:.1f}" r="4" fill="#2c3e50"/>')
        svg_parts.append(f'<text x="{col_x["drug"]}" y="{y - 8:.1f}" text-anchor="middle">{label}</text>')
    for label, y in regimen_y.items():
        svg_parts.append(f'<circle cx="{col_x["regimen"]}" cy="{y:.1f}" r="4" fill="#2c3e50"/>')
        short = label if len(label) < 40 else label[:37] + "..."
        svg_parts.append(f'<text x="{col_x["regimen"] + 10}" y="{y:.1f}" text-anchor="start" dominant-baseline="middle">{short}</text>')

    svg_parts.append("</svg>")
    svg = "\n".join(svg_parts)

    with open(output_path, "w") as f:
        f.write(svg)
    return output_path


# =====================================================================
# FULL PIPELINE ORCHESTRATION
# =====================================================================

def run_full_pipeline(
    genes: List[str],
    curated_gene_drugs: Dict[str, List[str]],
    real_gene_drugs: Dict[str, List[str]],
    patient_barcodes: List[str],
    maf_glob_pattern: str,
    gene_panel: List[str],
    resolve_tp53_fn: Callable,
    drug_graph_cls,
    synergy_net_cls,
    hcos_fn: Callable,
    mdcoe_fn: Callable,
    config: PipelineConfig,
    compute_diff: bool = True,
    export: bool = True,
    generate_visuals: bool = True,
) -> Dict:
    """Full loop: validate schemas -> parallel literature discovery ->
    provenance-tracked merge -> cohort-wide scoring -> frequency
    aggregation -> (optional) before/after diff -> (optional) export ->
    (optional) visual summary."""
    validate_gene_drug_schema(curated_gene_drugs, "curated_gene_drugs")
    validate_gene_drug_schema(real_gene_drugs, "real_gene_drugs")

    logger.info(f"pipeline start version={PIPELINE_VERSION} git={get_git_commit()}")

    discovery_result = parallel_discovery(
        genes, curated_gene_drugs, real_gene_drugs,
        max_workers=config.max_workers,
        cancer_context=config.cancer_context,
        max_papers_per_gene=config.max_papers_per_gene,
        api_key=config.api_key,
    )
    merged_gene_drugs, provenance = build_merged_pool_with_provenance(discovery_result, curated_gene_drugs)

    cohort_result = run_cohort_wide_analysis(
        patient_barcodes, maf_glob_pattern, gene_panel, merged_gene_drugs,
        resolve_tp53_fn, drug_graph_cls, synergy_net_cls, hcos_fn, mdcoe_fn,
        max_patients=config.max_patients,
    )
    freq = aggregate_regimen_frequency(cohort_result)

    diff = None
    if compute_diff:
        diff = diff_regimens_with_and_without_discovery(
            patient_barcodes, maf_glob_pattern, gene_panel,
            curated_gene_drugs, merged_gene_drugs,
            resolve_tp53_fn, drug_graph_cls, synergy_net_cls, hcos_fn, mdcoe_fn,
            max_patients=config.max_patients,
        )

    export_paths, visual_paths = {}, {}
    if export:
        export_paths = export_results(config.output_dir, discovery_result, merged_gene_drugs, provenance, cohort_result, freq, diff)
    if generate_visuals:
        os.makedirs(config.output_dir, exist_ok=True)
        visual_paths["bar_chart"] = generate_frequency_bar_chart(freq, os.path.join(config.output_dir, "regimen_frequency_bar.png"))
        visual_paths["flow_svg"] = generate_gene_drug_regimen_flow_svg(provenance, cohort_result, os.path.join(config.output_dir, "gene_drug_regimen_flow.svg"))

    logger.info("pipeline complete")
    return {
        "discovery_result": discovery_result,
        "merged_gene_drugs": merged_gene_drugs,
        "provenance": provenance,
        "cohort_result": cohort_result,
        "freq": freq,
        "diff": diff,
        "export_paths": export_paths,
        "visual_paths": visual_paths,
    }
