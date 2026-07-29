"""
Reproducibility metadata: pipeline version, git commit, and run timestamp
stamped onto every output artifact.
"""

from __future__ import annotations

import os
import subprocess
from datetime import datetime, timezone
from typing import Optional

import pandas as pd

PIPELINE_VERSION = "3.0.0"


def get_git_commit(cwd: Optional[str] = None) -> str:
    """Returns the short git commit hash of cwd (defaults to this file's
    directory), or a clearly-labeled fallback string if cwd isn't a git
    repo or git isn't available -- never raises, since a missing commit
    hash shouldn't crash a scientific run over a traceability nicety."""
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=cwd or os.path.dirname(os.path.abspath(__file__)),
            capture_output=True, text=True, timeout=5,
        )
        if result.returncode == 0:
            return result.stdout.strip()
        return "unknown (not a git repository)"
    except (FileNotFoundError, subprocess.SubprocessError):
        return "unknown (git unavailable)"


def stamp_dataframe(df: pd.DataFrame, cwd: Optional[str] = None) -> pd.DataFrame:
    """Returns a COPY of df with three added columns for traceability.
    Never mutates the input."""
    out = df.copy()
    out["_pipeline_version"] = PIPELINE_VERSION
    out["_git_commit"] = get_git_commit(cwd)
    out["_run_timestamp_utc"] = datetime.now(timezone.utc).isoformat()
    return out
