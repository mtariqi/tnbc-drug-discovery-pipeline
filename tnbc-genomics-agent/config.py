"""
config.py
Centralised configuration for the TNBC Genomics Agent pipeline.
Override any value by setting the corresponding environment variable.
"""

import os

# ── Anthropic API ─────────────────────────────────────────────────────────────
ANTHROPIC_API_KEY = os.getenv("ANTHROPIC_API_KEY", "")
AI_MODEL          = os.getenv("TNBC_AI_MODEL", "claude-sonnet-4-6")
AI_MAX_TOKENS     = int(os.getenv("TNBC_AI_MAX_TOKENS", "1000"))

# ── Pipeline defaults ─────────────────────────────────────────────────────────
DEFAULT_VCF_PATH    = os.getenv("TNBC_VCF_PATH",    "data/sample/patient_1.vcf")
DEFAULT_OUTPUT_DIR  = os.getenv("TNBC_OUTPUT_DIR",  "reports")
PIPELINE_VERSION    = "1.0.0"

# ── Variant filtering thresholds ─────────────────────────────────────────────
MIN_ALLELE_FREQ     = float(os.getenv("TNBC_MIN_AF",       "0.05"))   # 5% VAF floor
MIN_READ_DEPTH      = int(os.getenv("TNBC_MIN_DEPTH",      "20"))     # min coverage
MIN_SEVERITY_SCORE  = float(os.getenv("TNBC_MIN_SEVERITY", "0.0"))    # include all by default

# ── Redundancy scoring ────────────────────────────────────────────────────────
REDUNDANCY_SIGMOID_FACTOR = 0.8   # controls steepness of pathway redundancy curve
HIGH_BYPASS_THRESHOLD     = 0.6   # score at/above which bypass risk is HIGH
MODERATE_BYPASS_THRESHOLD = 0.3   # score at/above which bypass risk is MODERATE

# ── Reporting ─────────────────────────────────────────────────────────────────
MAX_COMBO_SUGGESTIONS = 10   # max combination therapy suggestions in report
MAX_BYPASS_TARGETS    = 20   # max bypass-risk rows in report
