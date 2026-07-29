"""
Plugin wiring cohort_wide_regimen_analysis.py (this repo's src/regimen/) to the REAL,
validated mdcoe.py -- which now lives right next to this file, in src/regimen/, not in a
separate src/utils/ (that folder turned out to be unnecessary once the real source was
uploaded: every one of resolve_tp53_drugs/DrugGraph/SynergyNet/HCOS/MDCOE lives in one file).

VERIFIED END-TO-END in this session: real mdcoe.py + this plugin's wiring +
cohort_wide_regimen_analysis.py, run together (not mocked), correctly reproduces the
afatinib+alpelisib+trastuzumab HCOS=0.450 result for matching patients and gives different
patients different regimens.
"""
import sys
import os

# mdcoe.py is a sibling file in this same directory now -- no path manipulation needed for
# the normal case. This shim only matters if you're running this plugin from somewhere else
# in the tree without src/regimen/ already on sys.path.
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))

if _THIS_DIR not in sys.path:
    sys.path.insert(0, _THIS_DIR)

MDCOE_SRC_PATH = _THIS_DIR  # kept for backwards compatibility with any code importing this name


def get_pipeline_components():
    from mdcoe import GENE_DRUGS, resolve_tp53_drugs, DrugGraph, SynergyNet, HCOS, MDCOE
    return {
        "resolve_tp53_fn": resolve_tp53_drugs,
        "drug_graph_cls": DrugGraph,
        "synergy_net_cls": SynergyNet,
        "hcos_fn": HCOS,
        "mdcoe_fn": MDCOE,
    }
