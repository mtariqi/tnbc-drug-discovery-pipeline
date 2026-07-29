"""
Package-wide logging configuration.

Call configure_logging() once, early (e.g. in a CLI entrypoint or at the
top of a notebook/script), to attach a single formatted handler to the
"tnbc_regimen_pipeline" logger. Every module below gets its own logger via
`logging.getLogger(__name__)` and relies on propagation to this one
handler -- so log lines are automatically prefixed with the actual
submodule that produced them (e.g. "tnbc_regimen_pipeline.discovery.parallelization")
without each module needing its own handler setup.

If configure_logging() is never called, the standard library's default
"no handlers found" behavior applies (a one-time warning to stderr) -- the
package deliberately does not call basicConfig() itself at import time,
since libraries configuring logging on import is an anti-pattern that
fights whatever logging setup the importing application already has.
"""

from __future__ import annotations

import logging

PACKAGE_LOGGER_NAME = "tnbc_regimen_pipeline"


def configure_logging(level: int = logging.INFO, stream=None) -> logging.Logger:
    """Attaches one formatted StreamHandler to the package root logger.
    Idempotent: calling it again just adjusts the level rather than
    stacking duplicate handlers."""
    logger = logging.getLogger(PACKAGE_LOGGER_NAME)
    logger.setLevel(level)

    if not logger.handlers:
        handler = logging.StreamHandler(stream)
        handler.setFormatter(
            logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s")
        )
        logger.addHandler(handler)
        logger.propagate = False

    return logger
