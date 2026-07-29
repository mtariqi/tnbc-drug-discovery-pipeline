"""
Configuration management.

PipelineConfig centralizes parameters that were previously hard-coded at
call sites (max_papers_per_gene, cancer_context, etc.) so a run is
reproducible from a committed file. config.yaml shipped alongside this
module is the package's default; load_default_config() reads it via the
actual installed package location, so it resolves correctly whether the
package is run from a source checkout or installed via pip.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Optional, Union

import yaml

_THIS_DIR = Path(__file__).parent
DEFAULT_CONFIG_PATH = _THIS_DIR / "config.yaml"


@dataclass
class PipelineConfig:
    cancer_context: str = "triple negative breast cancer"
    max_papers_per_gene: int = 10
    api_key: Optional[str] = None
    max_patients: Optional[int] = None
    max_workers: int = 8
    cache_dir: str = "./.pipeline_cache"
    output_dir: str = "./outputs"
    min_candidates_for_gap: int = 2


def load_config(path: Union[str, Path]) -> PipelineConfig:
    """Loads a PipelineConfig from a .yaml/.yml or .json file. Raises
    TypeError if the file contains a key PipelineConfig doesn't define --
    deliberately, so a misspelled config key fails loudly at load time
    rather than being silently ignored at run time."""
    path = str(path)
    with open(path) as f:
        if path.endswith((".yaml", ".yml")):
            data = yaml.safe_load(f) or {}
        elif path.endswith(".json"):
            data = json.load(f)
        else:
            raise ValueError(f"Unrecognized config extension for {path!r}; use .yaml, .yml, or .json")
    return PipelineConfig(**data)


def load_default_config() -> PipelineConfig:
    """Loads the package's shipped default config (config/config.yaml)."""
    return load_config(DEFAULT_CONFIG_PATH)


def save_config(config: PipelineConfig, path: Union[str, Path]) -> None:
    """Writes a PipelineConfig back out, so a run's exact parameters can be
    committed alongside its results for reproducibility."""
    path = str(path)
    data = asdict(config)
    with open(path, "w") as f:
        if path.endswith((".yaml", ".yml")):
            yaml.safe_dump(data, f, sort_keys=False)
        elif path.endswith(".json"):
            json.dump(data, f, indent=2)
        else:
            raise ValueError(f"Unrecognized config extension for {path!r}; use .yaml, .yml, or .json")
