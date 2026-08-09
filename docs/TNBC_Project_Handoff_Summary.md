# TNBC Manuscript Project — Handoff Summary
*For continuity in a new chat. Paste or reference this at the start of a new conversation, alongside re-uploading the real files listed below.*

## Current Deliverables (real files, already saved)
- **`TNBC_Manuscript_Draft.docx`** — the current, real manuscript (~26 pages, 7 figures, 33 references)
- **`TNBC_Professor_Deck.pptx`** — 18-slide walkthrough deck for faculty review
- **`TNBC_QA_Prep_Deck.pptx`** — 16-slide anticipated-questions defense prep deck
- **`extend_essentiality_with_sanger_score.py`** — real, tested script for Sanger Project Score cross-validation (independent-arm design, not pooling)
- **`run_real_essentiality_comparison.py`** — real, tested orchestration script that independently reproduced the manuscript's R²=0.337/−0.029 essentiality finding from raw DepMap files
- **`build_drug_repurposing_table.py`** — real script producing the drug-repurposing table from real DGIdb data

## Manuscript Structure (current)
1. Introduction (1.1 Related Work, 1.2 Contributions — 8 bullets)
2. Materials and Methods (2.1–2.10)
3. Experimental Setup (3.1–3.4, including the FunMap ablation design)
4. Results (4.1–4.9)
5. Discussion (5.1–5.7)
6. Limitations
7. Conclusions
References [1]–[33]

## The Core Scientific Framework
- **CTS** = 0.30·Centrality + 0.25·Essentiality + 0.25·Survival + 0.20·Druggability, over a real 90-gene RTK/NRTK panel
- **PairCTS**, **TripletCTS** — hierarchical extensions for pair/triplet prioritization (real equations in Methods §2.2–2.4)
- **HCOS/MDCOE** — a **separate**, patient-specific heuristic (confirmed via real repository-wide code search) that produced the focal-patient headline result (afatinib+alpelisib/capivasertib+trastuzumab, HCOS=0.450) — explicitly **not** built on CTS/PairCTS/TripletCTS, a major self-correction documented in §4.3/§5.3
- **Focal patient**: TCGA-AO-A128 (HER2-mutant, PTEN-null, TNBC)

## Major Findings, Verified Against Real Data
1. **DepMap essentiality validation**: gene identity alone (R²=0.337) beats a trained per-sample multi-omic model (R²=−0.029) — **independently re-derived from raw DepMap files this session** after the original computing script could not be located anywhere in the codebase despite exhaustive search. Real provenance: 34 TNBC lines identified by subtype label → 9 lack complete dependency data → 25 analysis-ready lines used for the real R² computation.
2. **CPTAC-STRING cross-validation**: 38/105 significant real correlations matched STRING edges; ERBB3↔ABL1 (r=0.506) is the strongest unmatched, independently literature-confirmed.
3. **FunMap vs. STRING comparison** (§4.7): 13% edge overlap; median rank shift 117 positions; only 5/20 top pairs shared despite high aggregate Spearman correlation (ρ≈0.93–0.98) — the "high correlation masks reordering" finding, now the manuscript's flagship methodological argument (Figure 7, "Evidence-Driven CTS Design").
4. **Track C** (TNBC molecular subtyping): tested as a candidate 5th CTS layer, excluded on real evidence (silhouette 0.070/0.058/0.025 across 3 clustering methods) — reframed as "a feature-selection result, not a failed experiment," arguably the strongest sentence in the manuscript.

## Major Self-Corrections Made This Project (all disclosed in the manuscript itself)
- **mdcoe/CTS attribution**: found via real repository-wide code search that the headline regimen result does NOT come from CTS/PairCTS/TripletCTS as originally implied — corrected in §4.3, §5.3, and reflected in Methods.
- **Redundancy-penalty over-application bug**: applying the top-N-only redundancy penalty across the *entire* ranking produced 15–35× inflated penalties; caught and fixed, now documented as a methodological lesson in §3.4.
- **Reference list errors**: 6 of the original 18 references had real errors (wrong author, wrong journal, wrong year, wrong first-author attribution) — all independently verified and corrected via real web search.
- **34→25 cell line provenance**: exhaustively traced why the manuscript says "25" TNBC lines when the source identification script produces 34 — resolved via real, from-scratch reproduction (see Major Findings #1).

## Professor/Reviewer Feedback Incorporated
Two real papers were reviewed and partially incorporated per multi-round professor feedback:
- Li et al. 2026 (*Discover Oncology*) — TNBC vs. non-TNBC mutational comparison; used to correct PAM50≠TNBC conflation (§2.6) and to add a carefully-caveated RTK-RAS citation (§1) — explicitly NOT claiming RTK-RAS is more altered in TNBC than non-TNBC (real numbers: 55.3% TNBC vs. 89.4% non-TNBC — the opposite direction).
- Chappell et al. 2021 (*Molecular Omics*) — TNBC multi-omics integration; used to qualify the PTEN-loss→PI3K/AKT inference in §4.3/§5.1 (PTEN-null HCC1937 showed *downregulated* PI3K/AKT/mTOR pathway activity) and to add a new genomic-state-vs-functional-state paragraph to the Introduction.
- **Title change ("multi-omic"→"multi-evidence") was explicitly considered and declined** — current title retained pending explicit sign-off from Dr. Miah/coauthors.

## Established Working Discipline (carry this forward)
- **Never fabricate a number.** Every statistic in the manuscript is either real (verified against actual data/code) or explicitly flagged as a documentation estimate needing verification.
- **Verify claims against real files/code before accepting them**, including your own or a professor's paraphrasing of "the manuscript currently says X."
- **Distinguish real computed results from citations of results** — the whole 34→25 investigation happened because a number was being cited without ever being traced to its actual computation.
- User's real environment: `~/projects/tnbc-genomics-agent-cloud/` and `~/rtk_nrtk_tnbc/`, fish shell, conda environments (mix-ups with wrong environment/missing packages have happened multiple times — check `which python3` / `conda env list` if errors occur).

## Genuinely Open Threads (not yet done)
- Real Sanger Project Score cross-validation (script ready, not yet run for real — see `extend_essentiality_with_sanger_score.py`)
- Kinase-redundancy paragraph from Chappell et al. (identified as "if space permits" tier, not yet added)
- Stronger CPTAC/FunMap interpretation paragraph (same tier)
- Paper 2 roadmap: FunMap-aware TripletCTS extension, real drug-response validation (GDSC/CTRP), full 93-patient cohort regimen ranking, phosphoproteomic CTS extension, mutation-context MDCOE features
- Final read-through by Dr. Miah still pending as of this summary

## How to Resume in a New Chat
1. Re-upload `TNBC_Manuscript_Draft.docx` (and any script you want to continue) so the new chat has the real, current file rather than a memory of it.
2. Paste or attach this summary.
3. State what you want to work on next.
