"""
Command-line entrypoint for tnbc_regimen_pipeline.

HONEST DESIGN NOTE ON WHAT A CLI CAN AND CAN'T DO HERE: the full pipeline
needs five callables/classes from your existing, separately-validated
MDCOE/HCOS scoring code (resolve_tp53_fn, DrugGraph, SynergyNet, hcos_fn,
mdcoe_fn) -- these are Python objects, not something a command line or a
YAML file can express directly. Rather than reimplementing or stubbing
them, every command here that needs them takes a --plugin argument (or a
'plugin:' field in a job-config file): a dotted "module:function" path to
a function in YOUR code that returns a dict with those five keys. This
keeps the actual scoring logic exactly where it's already validated,
outside this package, while still letting the whole pipeline run
end-to-end from a terminal.

Example plugin module (you write this, not part of this package):

    # my_mdcoe_plugin.py
    def get_pipeline_components():
        from my_validated_mdcoe_module import DrugGraph, SynergyNet, hcos, mdcoe, resolve_tp53
        return {
            "resolve_tp53_fn": resolve_tp53,
            "drug_graph_cls": DrugGraph,
            "synergy_net_cls": SynergyNet,
            "hcos_fn": hcos,
            "mdcoe_fn": mdcoe,
        }

    tnbc-pipeline full-pipeline --plugin my_mdcoe_plugin:get_pipeline_components ...
    tnbc-pipeline run job_config.yaml   # equivalent, single-file form; see run_cmd()
"""

from __future__ import annotations

import importlib
import json
import logging
import os
from typing import Callable, Dict

import click
import yaml

from ..config import load_config, load_default_config
from ..discovery.parallelization import parallel_discovery
from ..discovery.provenance import build_merged_pool_with_provenance
from ..pipeline.full_pipeline import run_full_pipeline
from ..utils.logging_setup import configure_logging

REQUIRED_PLUGIN_KEYS = {"resolve_tp53_fn", "drug_graph_cls", "synergy_net_cls", "hcos_fn", "mdcoe_fn"}


def _load_plugin(dotted_path: str) -> Dict:
    """Loads 'module.path:function_name', calls it, and validates the
    returned dict has exactly the keys the pipeline needs. Raises
    click.ClickException (not a raw traceback) on any failure, since this
    runs from a terminal."""
    if ":" not in dotted_path:
        raise click.ClickException(
            f"plugin must be in 'module.path:function_name' form, got {dotted_path!r}"
        )
    module_path, func_name = dotted_path.split(":", 1)
    try:
        module = importlib.import_module(module_path)
    except ImportError as e:
        raise click.ClickException(f"could not import plugin module {module_path!r}: {e}")
    try:
        factory: Callable = getattr(module, func_name)
    except AttributeError:
        raise click.ClickException(f"plugin module {module_path!r} has no attribute {func_name!r}")

    components = factory()
    missing = REQUIRED_PLUGIN_KEYS - set(components.keys())
    if missing:
        raise click.ClickException(
            f"plugin {dotted_path!r} returned a dict missing required keys: {sorted(missing)}"
        )
    return components


def _load_json_arg(path: str, name: str) -> Dict:
    try:
        with open(path) as f:
            return json.load(f)
    except (OSError, json.JSONDecodeError) as e:
        raise click.ClickException(f"could not read {name} from {path!r}: {e}")


def _resolve_inline_or_file(job: Dict, inline_key: str, file_key: str, name: str):
    """Job-config fields can be given inline (small dicts/lists directly
    in the YAML) or as a path to a separate JSON/text file (for large
    ones) -- not both."""
    if inline_key in job and file_key in job:
        raise click.ClickException(f"job config: specify either {inline_key!r} or {file_key!r}, not both")
    if inline_key in job:
        return job[inline_key]
    if file_key in job:
        path = job[file_key]
        if file_key.endswith("_file") and not path.endswith(".json"):
            with open(path) as f:
                return [line.strip() for line in f if line.strip()]
        return _load_json_arg(path, name)
    raise click.ClickException(f"job config must specify either {inline_key!r} or {file_key!r}")


def _load_job_config(path: str) -> Dict:
    """Parses a single 'run' job-config YAML. Schema:

        genes: [EGFR, PTEN, TYK2]
        gene_panel: [EGFR, PTEN, TYK2]
        curated_gene_drugs: {EGFR: [afatinib], ...}   # or curated_gene_drugs_json: path/to.json
        real_gene_drugs: {...}                          # or real_gene_drugs_json: path/to.json
        patient_barcodes: [PATIENT-001, ...]            # or patients_file: path/to.txt
        maf_glob_pattern: "data/*.maf.gz"
        plugin: my_mdcoe_plugin:get_pipeline_components  # dotted path -- see module docstring
        pipeline_config: path/to/pipeline_config.yaml    # optional, defaults to the packaged default
        output_dir: outputs/                             # optional override
        compute_diff: true                                # optional, default true
        generate_visuals: true                             # optional, default true

    This exists specifically because a pipeline "run" needs BOTH
    serializable data (genes, patients, file globs -- fine in YAML) AND
    non-serializable Python objects (DrugGraph, hcos_fn, etc. -- NOT
    representable in YAML). 'plugin' is how the latter gets in without
    pretending YAML can hold live code.
    """
    with open(path) as f:
        job = yaml.safe_load(f) or {}

    required = {"genes", "gene_panel", "maf_glob_pattern", "plugin"}
    missing = required - set(job.keys())
    if missing:
        raise click.ClickException(f"job config {path!r} missing required field(s): {sorted(missing)}")

    return job


@click.group()
@click.option("--verbose", is_flag=True, help="Enable INFO-level logging.")
def cli(verbose: bool):
    """tnbc-pipeline: TNBC combination-regimen discovery and scoring."""
    configure_logging(level=logging.INFO if verbose else logging.WARNING)


@cli.command("discover")
@click.option("--genes", required=True, help="Comma-separated gene list, e.g. EGFR,PTEN,TYK2")
@click.option("--curated-json", required=True, type=click.Path(exists=True), help="Path to curated_gene_drugs JSON")
@click.option("--real-json", required=True, type=click.Path(exists=True), help="Path to real_gene_drugs JSON")
@click.option("--config", "config_path", type=click.Path(exists=True), default=None, help="Path to a config YAML/JSON (defaults to the packaged default)")
@click.option("--output-dir", default=None, help="Overrides config.output_dir if given")
def discover(genes: str, curated_json: str, real_json: str, config_path: str, output_dir: str):
    """Run literature discovery + provenance-tracked merge only (no cohort scoring)."""
    config = load_config(config_path) if config_path else load_default_config()
    if output_dir:
        config.output_dir = output_dir

    curated_gene_drugs = _load_json_arg(curated_json, "curated_gene_drugs")
    real_gene_drugs = _load_json_arg(real_json, "real_gene_drugs")
    gene_list = [g.strip() for g in genes.split(",") if g.strip()]

    discovery_result = parallel_discovery(
        gene_list, curated_gene_drugs, real_gene_drugs,
        max_workers=config.max_workers,
        cancer_context=config.cancer_context,
        max_papers_per_gene=config.max_papers_per_gene,
        api_key=config.api_key,
    )
    merged, provenance = build_merged_pool_with_provenance(discovery_result, curated_gene_drugs)

    os.makedirs(config.output_dir, exist_ok=True)
    discovery_result.to_csv(os.path.join(config.output_dir, "discovery_result.csv"), index=False)
    with open(os.path.join(config.output_dir, "merged_gene_drugs.json"), "w") as f:
        json.dump(merged, f, indent=2)
    with open(os.path.join(config.output_dir, "provenance.json"), "w") as f:
        json.dump(provenance, f, indent=2)

    click.echo(f"Discovery complete: {len(discovery_result)} candidate extraction(s), "
               f"results written to {config.output_dir}")


@cli.command("full-pipeline")
@click.option("--genes", required=True, help="Comma-separated gene list")
@click.option("--curated-json", required=True, type=click.Path(exists=True))
@click.option("--real-json", required=True, type=click.Path(exists=True))
@click.option("--patients-file", required=True, type=click.Path(exists=True), help="One patient barcode per line")
@click.option("--maf-glob", required=True, help="Glob pattern for gzipped MAF files")
@click.option("--gene-panel", required=True, help="Comma-separated gene panel")
@click.option("--plugin", required=True, help="dotted 'module:function' returning the MDCOE/HCOS component dict (see module docstring)")
@click.option("--config", "config_path", type=click.Path(exists=True), default=None)
@click.option("--output-dir", default=None)
@click.option("--no-diff", is_flag=True, help="Skip the before/after discovery diff")
@click.option("--no-visuals", is_flag=True, help="Skip bar chart / flow SVG generation")
def full_pipeline_cmd(genes, curated_json, real_json, patients_file, maf_glob, gene_panel, plugin, config_path, output_dir, no_diff, no_visuals):
    """Run the full discovery -> cohort scoring -> diff -> export -> visuals pipeline, from individual flags."""
    config = load_config(config_path) if config_path else load_default_config()
    if output_dir:
        config.output_dir = output_dir

    curated_gene_drugs = _load_json_arg(curated_json, "curated_gene_drugs")
    real_gene_drugs = _load_json_arg(real_json, "real_gene_drugs")
    gene_list = [g.strip() for g in genes.split(",") if g.strip()]
    gene_panel_list = [g.strip() for g in gene_panel.split(",") if g.strip()]
    with open(patients_file) as f:
        patient_barcodes = [line.strip() for line in f if line.strip()]

    components = _load_plugin(plugin)

    result = run_full_pipeline(
        genes=gene_list,
        curated_gene_drugs=curated_gene_drugs,
        real_gene_drugs=real_gene_drugs,
        patient_barcodes=patient_barcodes,
        maf_glob_pattern=maf_glob,
        gene_panel=gene_panel_list,
        resolve_tp53_fn=components["resolve_tp53_fn"],
        drug_graph_cls=components["drug_graph_cls"],
        synergy_net_cls=components["synergy_net_cls"],
        hcos_fn=components["hcos_fn"],
        mdcoe_fn=components["mdcoe_fn"],
        config=config,
        compute_diff=not no_diff,
        export=True,
        generate_visuals=not no_visuals,
    )

    n_changed = 0 if result["diff"] is None else len(result["diff"])
    click.echo(
        f"Pipeline complete: {len(patient_barcodes)} patient(s) scored, "
        f"{n_changed} patient(s) changed by discovery, "
        f"results written to {config.output_dir}"
    )


@cli.command("run")
@click.argument("job_config_path", type=click.Path(exists=True))
def run_cmd(job_config_path: str):
    """Run the full pipeline from a single YAML job-config file (genes,
    patients, MAF glob, and a plugin dotted-path for the MDCOE/HCOS
    components -- see _load_job_config()'s docstring for the schema)."""
    job = _load_job_config(job_config_path)

    curated_gene_drugs = _resolve_inline_or_file(job, "curated_gene_drugs", "curated_gene_drugs_json", "curated_gene_drugs")
    real_gene_drugs = _resolve_inline_or_file(job, "real_gene_drugs", "real_gene_drugs_json", "real_gene_drugs")
    patient_barcodes = _resolve_inline_or_file(job, "patient_barcodes", "patients_file", "patient_barcodes")

    pipeline_config_path = job.get("pipeline_config")
    config = load_config(pipeline_config_path) if pipeline_config_path else load_default_config()
    if job.get("output_dir"):
        config.output_dir = job["output_dir"]

    components = _load_plugin(job["plugin"])

    result = run_full_pipeline(
        genes=job["genes"],
        curated_gene_drugs=curated_gene_drugs,
        real_gene_drugs=real_gene_drugs,
        patient_barcodes=patient_barcodes,
        maf_glob_pattern=job["maf_glob_pattern"],
        gene_panel=job["gene_panel"],
        resolve_tp53_fn=components["resolve_tp53_fn"],
        drug_graph_cls=components["drug_graph_cls"],
        synergy_net_cls=components["synergy_net_cls"],
        hcos_fn=components["hcos_fn"],
        mdcoe_fn=components["mdcoe_fn"],
        config=config,
        compute_diff=job.get("compute_diff", True),
        export=True,
        generate_visuals=job.get("generate_visuals", True),
    )

    click.echo("=== Discovery ===")
    click.echo(result["discovery_result"].to_string(index=False) if not result["discovery_result"].empty else "(no coverage gaps found)")
    click.echo("\n=== Merged gene->drug pool ===")
    click.echo(json.dumps(result["merged_gene_drugs"], indent=2))
    click.echo("\n=== Cohort ===")
    click.echo(result["cohort_result"].to_string(index=False))
    click.echo("\n=== Frequency ===")
    click.echo(result["freq"].to_string(index=False))
    if result["diff"] is not None:
        click.echo("\n=== Diff (patients changed by discovery) ===")
        click.echo(result["diff"].to_string(index=False) if not result["diff"].empty else "(no patient's top regimen changed)")
    click.echo(f"\nResults written to {config.output_dir}")


def main():
    cli()


if __name__ == "__main__":
    main()
