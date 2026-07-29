"""
tnbc_regimen_pipeline
=====================

Kinase-alteration-driven combination therapy scoring for TNBC (TCGA-BRCA
cohort), extended with literature-mining discovery of novel drug
candidates for under-covered genes.

Subpackages:
    discovery -- PubMed search, INN-suffix extraction, mandatory DGIdb
                 confirmation, caching, parallelization, provenance
    cohort    -- TNBC patient identification, MAF reading, cohort-wide
                 MDCOE/HCOS scoring, frequency aggregation
    pipeline  -- orchestration wiring discovery + cohort together, diffing,
                 exports, visual summary, reproducibility stamping
    config    -- PipelineConfig and YAML/JSON loaders
    cli       -- command-line entrypoint
    utils     -- logging setup, schema validation
"""

from .pipeline.reproducibility import PIPELINE_VERSION as __version__

__all__ = ["__version__"]
