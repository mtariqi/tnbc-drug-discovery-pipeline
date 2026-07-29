# TNBC Drug Regimen Discovery

![Banner](assets/banner.png)

[![Python 3.9+](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Search Verified](https://img.shields.io/badge/beam%20search-bug%20fixed%20%26%20re--verified-brightgreen)](docs/Ranking_Methodology_Report.docx)
[![Result](https://img.shields.io/badge/top%20HCOS-0.450-informational)](results/top10_ranked_regimens_TCGA-AO-A128.csv)
[![Status](https://img.shields.io/badge/status-active-success)]()

A beam-search + heuristic-scoring engine (**MDCOE** / **HCOS**) that ranks candidate multi-drug regimens against a real patient's genomic profile — including a real algorithmic bug found, fixed, and re-verified during this session.

---

## Key Result

For confirmed-TNBC focal patient **TCGA-AO-A128** (6 altered druggable genes: EGFR, PTEN, FLT1, TP53, ERBB2, TYK2 → 15 candidate drugs), searched across all regimen sizes (2–5 drugs):

| Rank | HCOS | Size | Regimen |
|---|---|---|---|
| 1 (tie) | **0.450** | 3-drug | **afatinib + alpelisib + trastuzumab** |
| 1 (tie) | **0.450** | 3-drug | afatinib + capivasertib + trastuzumab |
| 3–10 | 0.400 | 2-drug | (see [full results](results/top10_ranked_regimens_TCGA-AO-A128.csv)) |

**Why the top regimen makes sense:** it covers 3 of the patient's 6 altered genes — EGFR (afatinib), PTEN-loss-driven PI3K reactivation (alpelisib), and an ERBB2 activating mutation (trastuzumab) — all within one coherent signalling module, with real clinical precedent (the Phase II SUMMIT trial tested HER2-directed therapy in exactly this non-amplified, HER2-mutant profile).

**Why the exact tie matters:** alpelisib and capivasertib both act on the same PI3K/AKT consequence of this patient's PTEN loss. The scoring system reasons at the *pathway* level, not the individual-drug level — so it can't distinguish between two drugs hitting the same node. That's a real, informative limitation, not noise.

Full methodology, worked example, and component-level breakdown: [`docs/Ranking_Methodology_Report.docx`](docs/Ranking_Methodology_Report.docx)

---

## The Bug: Beam Search Only Returned One Fixed Regimen Size

The search algorithm was found to only return regimens of *exactly* one size — whichever maximum depth was specified — because it kept extending every surviving candidate forward at each step and only examined the final depth's survivors at the end. A genuinely optimal 3-drug regimen could be generated partway through the search, then extended into a lower-scoring 4- or 5-drug version, and never reported at all.

**Fix:** every candidate regimen generated at *any* depth (size 2 and up) is now retained in a single pool, and the final ranking is taken across the whole pool — not just the last depth explored.

**Verification, not assumption:** re-ran the fixed search at `max_depth=5` and confirmed it finds the known-correct 3-drug answer (HCOS = 0.450) automatically, without needing to already know the right depth to search at. Four additional regression checks confirmed the fix didn't break deduplication, result sorting, minimum-size filtering, or `top_k` limiting.

```python
from src.mdcoe import GENE_DRUGS, resolve_tp53_drugs, DrugGraph, SynergyNet, HCOS, MDCOE

genes = ["EGFR", "PTEN", "FLT1", "TP53", "ERBB2", "TYK2"]
tp53_drugs = resolve_tp53_drugs("Nonsense_Mutation")  # confirm against real MAF classification before trusting this
drugs = sorted(set(GENE_DRUGS.get(g, []) if g != "TP53" else tp53_drugs for g in genes))

graph = DrugGraph(drugs)
net = SynergyNet()
results = MDCOE(graph, net, HCOS, beam_width=50, max_depth=5, top_k=10)
```

---

## Scoring Design (HCOS)

```
HCOS = synergy + evidence − toxicity_penalty − diversity_penalty + size_bonus
```

- **Synergy** — pathway overlap between target genes, an RTK/intracellular complementarity bonus, plus a small literature-backed synergy table for specific known pairs.
- **Evidence** — extra credit when a regimen contains a pair with direct literature/trial precedent.
- **Toxicity penalty** — subtracts for shared toxicity categories across the regimen's drugs (sourced from real drug-label data).
- **Diversity penalty** — subtracts when multiple drugs in a regimen redundantly re-target the same gene.
- **Size bonus** — small, deliberately mild reward for larger regimens — mild enough that, as the results above show, several 2-drug regimens still outscore larger alternatives when the extra drug doesn't add real value.

---

## Important Caveats

- **HCOS is a documented placeholder heuristic**, not a trained or clinically validated model — explicitly stated in the source code itself. It reasons over pathway membership and a small hardcoded synergy table (currently 4 literature-backed pairs), not real experimental synergy screens.
- **The TP53 alteration classification** used in the example above (`Nonsense_Mutation` → nonsense-mediated decay) is carried over from an illustrative example in this project and has not yet been confirmed against this patient's actual MAF-reported classification. **(Note: this was subsequently confirmed — see `results/reports/TP53_Confirmation_Addendum.docx` in the consolidated repo.)**
- **This is a computational prioritization for wet-lab testing**, not an experimentally or clinically validated result.

---

## Repository Structure

```
.
├── src/
│   └── mdcoe.py                                    # Bug-fixed beam search + HCOS scoring
├── results/
│   └── top10_ranked_regimens_TCGA-AO-A128.csv      # Validated output
├── docs/
│   └── Ranking_Methodology_Report.docx             # Full methodology + worked example
├── requirements.txt
└── README.md
```

## Installation & Usage

```bash
pip install -r requirements.txt
python3 src/mdcoe.py   # runs the built-in demo cases
```

See `src/mdcoe.py`'s module docstring for the full API (`GENE_DRUGS`, `resolve_tp53_drugs`, `DrugGraph`, `SynergyNet`, `HCOS`, `MDCOE`).

## Related Work

This regimen-ranking system (HCOS) is separate from, and not yet integrated with, a parallel STRING/DepMap/TCGA-driven kinase-scoring system (CTS/TripletCTS) developed in a companion project — see [tnbc-kinase-scoring-pipeline](../tnbc-kinase-scoring-pipeline).

**Repo-consolidation note:** `tnbc-kinase-scoring-pipeline`'s own README lists a *second* copy of `mdcoe.py` in its own `src/`. These are two copies of the same file, potentially at different bug-fix states — this repo's copy is the one documented above as having the beam-search depth bug found and fixed; confirm the kinase-scoring-pipeline copy matches before treating them as interchangeable. See `tnbc-genomics-agent/INTEGRATION_NOTES.md` item 2 in the consolidated repo for the same flag raised from the other direction.

## License

[MIT](LICENSE)
