"""
Local response cache for PubMed/DGIdb calls.

Wraps agentic_regimen_discovery's two network-touching functions
(search_pubmed_for_combination_therapy, confirm_candidate_via_dgidb) with
joblib.Memory disk caching, keyed on their arguments. Applied by
REPLACING the names in agentic_regimen_discovery's own module namespace
at call time -- the sibling module's source is never edited.

Why wrap rather than edit the source: a caching concern (repeat-run
performance) is unrelated to the correctness of the underlying network
calls, which is independently tested. Wrapping keeps that boundary clean
and makes the cache trivially removable (uninstall_cache()).
"""

from __future__ import annotations

from typing import Optional

from joblib import Memory

from . import agentic_regimen_discovery as ard

_original_search = ard.search_pubmed_for_combination_therapy
_original_confirm = ard.confirm_candidate_via_dgidb
_memory: Optional[Memory] = None


def install_cache(cache_dir: str = "./.pipeline_cache") -> Memory:
    """Replaces agentic_regimen_discovery's network functions with
    joblib-cached versions. Safe to call more than once (re-installing
    just repoints to a fresh Memory at the given dir). Returns the Memory
    object so callers can inspect/clear it (memory.clear())."""
    global _memory
    _memory = Memory(location=cache_dir, verbose=0)
    ard.search_pubmed_for_combination_therapy = _memory.cache(_original_search)
    ard.confirm_candidate_via_dgidb = _memory.cache(_original_confirm)
    return _memory


def uninstall_cache() -> None:
    """Restores the original, uncached functions."""
    ard.search_pubmed_for_combination_therapy = _original_search
    ard.confirm_candidate_via_dgidb = _original_confirm
