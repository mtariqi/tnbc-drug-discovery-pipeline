# Known Limitations and Open Gaps

Consolidated from `results/reports/TNBC_Consolidated_Report.docx` §6 and
`Session2_Report_DepMap_Extension.docx` §7, plus the new Track B gap identified this session.

## Treatment-response prediction (Track B) — not attempted, by design

The downloaded TCGA-BRCA clinical file has survival/vital-status fields only — confirmed
directly, not assumed. No pCR, RECIST, or drug-response field exists to train against.
Building this track from the current data would mean inventing a proxy label, which breaks
this project's real-data-only discipline (see the CTS/HCOS-gap handling below for how this
project treats a real gap: document it, don't paper over it).

**Real path forward:** a supplementary, TNBC-enriched cohort with documented neoadjuvant
response — GSE25066, GSE25055/GSE25065, or I-SPY2 (GSE194040) — as a separate data-acquisition
step, not a same-file shortcut.

## No real combination-synergy data

No pairwise drug-synergy dataset (GDSC2, PRISM) is integrated yet. PairCTS's synergy term and
TripletCTS's escape-route-closure term currently contribute no signal to any score.

## HCOS mechanism-gap

HCOS reasons at the kinase-pathway-overlap level and cannot yet credit a real,
literature-confirmed non-kinase-inhibitor synergy (demonstrated by the
carboplatin+eprenetapopt result ranking below the top kinase-inhibitor regimen despite real
supporting evidence — see `TNBC_Regimen_Recommendation_Final.docx` §4).

## Small ML training cohort (DepMap dependency model)

Trained on 25 confirmed-TNBC cell lines — enough to detect the reported gene-identity effect,
not powered to rule out a modest real per-sample multi-omic effect a larger cohort might
detect. No formal power calculation has been run to state a specific, defensible bound.

## Discovery-pipeline recall ceiling

The agentic literature-mining pipeline only recognizes standard-INN-named drugs already
present in DGIdb; it cannot surface a genuinely novel, unregistered compound.

## Two copies of mdcoe.py, possibly at different bug-fix states

Confirmed from the actual upstream READMEs (`docs/source_repo_readmes/`): both
`TNBC-drug-regimen-discovery` and `tnbc-kinase-scoring-pipeline` ship their own `mdcoe.py`.
The drug-regimen repo documents a real beam-search depth bug found and fixed there; the
kinase-scoring repo's own Known Limitations section separately states its copy "is not yet
reconciled" against CTS/TripletCTS. Not yet resolved — `src/regimen/my_mdcoe_plugin.py`
currently points at the drug-regimen-discovery copy specifically because of its documented fix,
but the two files have not been diffed against each other.

## CTS/HCOS and mdcoe.py's HCOS are two separate, uncross-validated scoring systems

Stated directly in `tnbc-kinase-scoring-pipeline`'s own README: "`mdcoe.py`'s HCOS scoring and
this repository's CTS/TripletCTS scoring are two separate systems for a related goal; they are
not yet reconciled or cross-validated against each other," and "`PairCTS`/`TripletCTS` have not
yet been run against this session's real, validated 90-kinase CTS output." Both remain open.

## gdc_local_data_loader.py's smoke test can't currently run standalone

Its own docstring states it "was tested against a synthetic directory tree built with the real
GDC STAR-counts TSV header format and real MAF column layout" — but the code that actually
*builds* that synthetic tree isn't present in the file (`_run_smoke_tests()` expects a
pre-existing `test_root` directory with a `gdc_sample_sheet*.tsv` and expression files already
in it, and crashes with `StopIteration` if that directory doesn't exist). Confirmed by actually
running it, not assumed from reading it. Rather than reconstruct a fake GDC fixture from
general knowledge of the STAR-Counts format (risky to guess at for a real genomics pipeline —
exactly the kind of silent assumption this project avoids elsewhere), follow the module's own
documented fix: run `inspect_gdc_directory("~/rtk_nrtk_tnbc/data/raw/tcga_brca")` on the real
data first and compare its output against what `_run_smoke_tests()` expects, or locate whatever
script originally generated that fixture directory during the real session.

## PairCTS/DrugPairScore ceiling effect from max-aggregation across targets — confirmed on real data

`best_target_pair_cts()` scores a drug pair by taking the **maximum** `pair_cts()` across every
combination of the two drugs' targets, not an average or profile-aware comparison. Run for
real against `~/rtk_nrtk_tnbc/data/`, this produces heavy exact ties at the top of the ranking
(`DrugPairScore = 0.656483` for most top-20 pairs) — confirmed by directly inspecting which
target-kinase-pair actually won the `max()` for several tied pairs: two independent drug pairs
(`dacomitinib`+`hesperadin` and `crizotinib`+`hesperadin`) both resolved to the literal same
winning kinase pair (`ERBB2`, `EPHA5`), same score, out of 273 and 507 possible combinations
respectively. `ERBB2` is the globally highest-CTS kinase in the 90-kinase panel — so any drug
pair where one drug's target set happens to include `ERBB2` and the other includes a
similarly-favorable kinase collapses onto the same ceiling score, regardless of what either
drug's broader target profile looks like. `hesperadin` alone contributing 273–507 target-pair
combinations (evidently a broad, non-selective research compound) makes this collision likely
for almost any pairing involving it.

**Second attempt: redundancy penalty for repeated winning kinase pairs — mixed, real result.**
`src/scoring/redundancy_penalized_pair_scoring.py` tracks which specific kinase pair "wins" the
max-search for each drug pair, then penalizes scores in proportion to how often that same
winning pair recurs elsewhere in the ranking (unique winning pairs get zero penalty). Tested
against the real 116-drug dataset at two penalty weights, with a genuine, honestly-mixed result:

| Metric | Original | weight=0.05 | weight=0.15 |
|---|---|---|---|
| Top-20 exact-score ties | 20/20 | 9/20 | **7/20** |
| Distinct winning kinase pairs in top 20 | ~1 | 12 | **13** |
| Unique scores, full 1770-pair ranking | 430 | 116 | **75** |

**Top-of-ranking diversity genuinely improved** (fewer ties, far more distinct kinase pairs
represented, monotonically better at higher weight) — `hesperadin` no longer dominates every
top result. **Full-ranking uniqueness genuinely got worse**, monotonically, at both tested
weights — not a tuning artifact fixable by picking a different weight, a real structural
trade-off. Plausible explanation (stated as a hypothesis, not verified with an actual score
histogram): `pair_cts()`'s raw scores are already fairly coarse/discretized (only 430 unique
values across 1770 pairs to start with), and summing a second, independently-discretized
quantity (integer frequency × fixed weight) on top can produce a coarser combined result in
the tail of the ranking even while separating the top.

**Decision made: adopted as the default for top-N reporting (weight=0.15), original kept
available for full-ranking use.** `run_pairs_and_triplets.py`'s printed "Top N pairs" and its
triplet-building drug pool now use `rank_all_pairs_redundancy_penalized()` by default (patched
via `apply_redundancy_penalty_default.py`). The plain, unpenalized `rank_all_pairs()` remains
importable and unchanged in `kinase_scoring_pipeline.py` for anyone who explicitly wants the
full 1770-pair ranking (e.g. via `compare_redundancy_penalty.py`) — nothing about the
underlying scoring formula changed, only which ranking `run_pairs_and_triplets.py` reports by
default.

**Follow-up, checked against real sources rather than assumed**: `altiratinib` and `metatinib`
appearing together across multiple top pairs/triplets (both converging on the same
`(EPHA5, KDR)`/`(KDR, MET)`-type winning kinase pair) was flagged as a possible duplicate-drug
or data-quality artifact, given how consistently they co-occurred. Verified via web search
(PMC8342607, PubMed 33749934): both are real, independently-developed compounds —
metatinib tromethamine is a distinct small-molecule c-MET/VEGFR2 inhibitor with its own
published Phase I trial, not a synonym or corrupted variant of altiratinib (a separate,
Deciphera-developed MET/TIE2/VEGFR2/TRKA inhibitor). Their consistent co-occurrence is real
signal, not a bug: two mechanistically similar real drugs correctly converging on the same
target pair. Worth keeping in mind when interpreting the ceiling effect above — some (not all)
of the observed convergence reflects genuine pharmacological similarity between distinct real
compounds, not purely an artifact of the max-aggregation scoring formula.

## Confirmed: compute_cts() genuinely reproducible from raw data, with one resolved methodological question

Built `build_and_verify_real_cts.py` to assemble the real `kinase_df` `compute_cts()` expects
directly from raw files (`kinase_90_list.txt`, DepMap essentiality, TCGA-BRCA survival, STRING
edges via the already-validated `compute_centrality()`, and DGIdb interactions) — rather than
only ever loading the already-computed `CTS_all_90_kinases.tsv` as `run_pairs_and_triplets.py`
does.

**Open question when first built:** how should raw DGIdb interaction rows aggregate into the
per-kinase `dgidb_score` druggability input? Two hypotheses were tested empirically against the
real, existing `CTS_all_90_kinases.tsv` as ground truth, rather than picking one and moving on:

| Aggregation | Mean abs diff (90 kinases) | Max abs diff | Kinases closer |
|---|---|---|---|
| `mean(interaction_score)` per kinase | 0.0548 | 0.166 | 21/90 |
| `count(distinct non-junk drugs)` per kinase | **0.0213** | **0.034** | **69/90** |

**Count of distinct drugs is the confirmed correct (or much closer) methodology** — a ~2.6x
reduction in total discrepancy, and the remaining max difference (0.034) is small enough to
plausibly be fully explained by the one still-unavailable input (`chembl_count`/`trial_stage`
data, which the original computation may have had access to and this reconstruction doesn't)
rather than a further methodological mismatch. `run_pairs_and_triplets.py`'s use of the
already-computed file remains valid and doesn't need to change; this was specifically about
being able to regenerate CTS from scratch if the underlying raw data ever changes.

## Corrected: no actual mdcoe.py duplication

An earlier version of this doc flagged two possibly-diverged copies of `mdcoe.py` across the
two original repos. **Checked against the real uploaded source — this doesn't hold.**
`tnbc-kinase-scoring-pipeline`'s real `src/` has no `mdcoe.py` at all; its README's listed
structure was stale. Only one real copy exists, in `TNBC-drug-regimen-discovery`, and it has
been verified end-to-end (see `src/regimen/` — the real `mdcoe.py` reproduces the
afatinib+alpelisib+trastuzumab, HCOS=0.450 result exactly, and wires correctly into
`cohort_wide_regimen_analysis.py` via `my_mdcoe_plugin.py`).

## R script (assemble_tcga_brca_data.R) crashes on real data -- GDCprepare() itself, not fixable downstream

Running the R assembly script for real (needed to get expression data for Track C) crashed
twice, consistently at the same point ("100% then crash") -- first attempt suspected the
crash was in the post-`GDCprepare()` transpose/data-frame conversion and added gene subsetting
before that step; this did NOT help, confirming the crash originates INSIDE `GDCprepare()`
itself while assembling the full ~1098-patient x ~60,000-gene `SummarizedExperiment` object,
before any downstream code (including the subsetting fix) ever runs. `dmesg` was blocked by
permissions so this couldn't be confirmed via kernel OOM logs directly, but the machine has
0B swap configured (`free -h`), which is consistent with abrupt OOM kills rather than graceful
slowdown.

**Real fix, reusing an already-built, already-documented solution to this exact problem**:
your own `tnbc-kinase-scoring-pipeline` README documents fixing this identical `GDCprepare()`
OOM crash once already, for the kinase panel specifically ("replaced with a file-by-file Python
loader that never holds more than one sample in memory"). `build_real_expression_for_track_c.py`
applies that same already-validated approach (`gdc_local_data_loader.py`'s
`build_expression_matrix()`, reading each raw STAR-counts TSV file directly, restricted
immediately to the ~126 needed genes) to Track C's marker genes, bypassing R and
`GDCprepare()` entirely rather than trying to make the R route survive.

**Also found via this real-data attempt**: `prepare_tnbc_expression_matrix()` had a real gap --
its zero-variance check (`var() == 0`) silently failed to catch an entirely-missing (all-NaN)
gene column, since `.var()` on all-NaN data returns `NaN`, not `0`. Fixed to explicitly drop any
gene with even one missing value (conservative -- not imputed) before NMF, which can't accept
NaN input at all.

## Track C real-cohort run: silhouette score much lower than synthetic testing -- diagnosed, two real fixes confirmed

Running Track C against the real 1095-patient/119-gene TNBC expression data gave a
silhouette score of only 0.070 (vs. 0.70+ on marker-only synthetic data) -- essentially no
real cluster separation. Two hypotheses were proposed and each CONFIRMED as a real, measurable
mechanism via controlled synthetic tests (not just asserted):

1. **Noise-gene dilution**: the real expression file combines 90 kinase-panel genes (reused
   from the CTS/HCOS scoring project) with only 29 additional genes actually relevant to
   Lehmann subtyping. Confirmed on synthetic data: silhouette 0.477 (90 noise genes included) →
   0.718 (restricted to marker genes only) -- a real, substantial improvement, not marginal.
2. **No log-transform**: raw TPM values are highly right-skewed; NMF run on raw values can be
   dominated by scale rather than genuine pattern. Confirmed: cluster purity 0.678 (unlogged,
   scale-skewed genes present) → 1.000 (log1p-transformed) -- a dramatic improvement.

**Both fixed and made the new defaults**: `prepare_tnbc_expression_matrix()` now takes
`restrict_to_genes` and `log_transform` parameters; `run_tnbc_subtyping_pipeline()` defaults to
`restrict_to_marker_genes=True` and `log_transform=True`. Also fixed in the same pass: a real
gap where an entirely-missing gene column (NaN variance, not 0) silently evaded the
zero-variance check and crashed NMF -- now explicitly dropped with a printed note.

**Re-run against the real cohort with both fixes: no real improvement.** Silhouette went
0.070 → 0.058 -- essentially unchanged, arguably slightly worse. This is a genuinely humbling
result: both mechanisms were real and confirmed on synthetic data, but that only establishes
they're real mechanisms, not that they're *sufficient* to fix real data -- the synthetic tests
used artificially large, clean effect sizes (a flat +8.0 boost) that don't represent how subtle
and noisy real biological signal actually is. Confirming a mechanism in a controlled toy case
doesn't guarantee it resolves the problem in practice, and here it didn't.

**Real path forward, built and confirmed on synthetic data**: rather than continue patching
the small curated marker panel, `gene_selection_mode="highly_variable"` clusters on the top
N most-variable genes across the FULL transcriptome (standard practice in real expression
clustering pipelines) instead of ~30 curated genes -- while marker-based characterization
still uses the fixed, curated marker sets afterward for interpretability, independent of what
genes clustering itself used. Confirmed on synthetic data (500-gene background, true signal
in only 13 genes): highly-variable selection correctly found 13/13 true informative genes,
clustering recovered 100% purity, and characterization correctly identified both clusters.

Getting the full transcriptome required a second real fix: `build_real_expression_for_track_c.py`
now supports a `--full` flag to fetch every gene (not just the 90-kinase+marker panel) via the
same memory-safe per-file loader -- confirmed this is NOT the same risk as `GDCprepare()`'s
crash, since a full ~1095 x ~20,000 matrix is only ~170MB as float64, nothing like what an
OOM-prone full `SummarizedExperiment` object required.

**Third attempt, full-transcriptome highly-variable genes: also weak, and worse, not better.**
Built the full transcriptome (1095 patients x 59,427 genes) via the memory-safe loader with a
`--full` flag (confirmed working at full scale: completed in under 5 minutes, no crash --
the real fix to the GDCprepare() problem holds up beyond the small gene-panel case it was
first proven on). Selected the top 2000 most-variable genes (standard practice) and
re-clustered: **silhouette = 0.025** -- worse than both marker-only attempts (0.070, 0.058),
not better.

**Synthesis across all three real attempts, and the actual conclusion:** three genuinely
different methodologies -- small curated marker panel (0.070), the same panel with
log-transform and noise-gene exclusion (0.058), and full-transcriptome highly-variable-gene
selection (0.025) -- all land in the same weak 0.02-0.07 silhouette range. When substantially
different approaches all converge on weak separation, that stops looking like a fixable
engineering problem specific to any one method and starts looking like a real property of the
cohort itself, not a bug to keep patching.

This is consistent with genuine, published findings about TNBC subtyping robustness, not just
a shortcoming of this pipeline: Lehmann's own follow-up work reduced the original 6 subtypes to
4 partly over reproducibility concerns, and subsequent literature has repeatedly noted TNBC
molecular subtypes show substantial continuous overlap rather than clean separability --
especially in bulk (not single-cell) RNA-seq data, which mixes tumor, stromal, and immune
signal together rather than isolating pure tumor expression. A silhouette score in the
0.02-0.15 range is not unusual in the real literature for this class of clustering problem.

**Recommendation: accept this as the honest, defensible result rather than continuing to
chase a higher number.** The low silhouette score most likely isn't hiding a fixable bug --
it may be accurately reflecting that this cohort's TNBC subtypes are genuinely not cleanly
discrete via bulk-expression clustering. The cluster/marker-characterization output
(`tnbc_subtype_assignments.tsv`) is still a legitimate result to report, just with this
honestly-earned caveat attached: treat the resulting subtype calls as soft, overlapping
tendencies, not confident discrete labels.

## CPTAC proteomics validation added, but not TNBC-restricted -- checked, not assumed

`cptac-proteomics-validation/` folds in real, measured CPTAC BRCA proteomics as an independent
check against STRING-*predicted* kinase crosstalk used in `PairCTS`/`TripletCTS`. Extensively
checked whether this cohort could be restricted to TNBC specifically before including it:
no ER/PR/HER2/PAM50 field exists anywhere in the `cptac` Python package's clinical data for
this cohort — checked `get_clinical(source='mssm')` (the only available clinical source for
BRCA), `get_follow-up()`, `get_medical_history()`, and the raw pan-cancer clinical file
directly. Real TNBC-restriction would require the original publication's supplementary tables,
not pursued further this session. This data is BRCA-wide (151 samples, all subtypes), not
TNBC-specific — every result from it should be read that way.

One real, precisely-scoped bug found and confirmed by actually re-running the pipeline (not
just reading the code): `correlation_analysis_final.py`'s plain-text summary prints raw,
truncated Ensembl compound IDs instead of gene symbols for its "strongest correlation" lines.
Confirmed cosmetic only — the underlying CSVs and the companion script
(`rtk_nrtk_network_plotly.py`) both correctly resolve real gene symbols, and both independently
agree on the same top result (`ERBB3 ↔ ABL1`, r=0.506).

**Follow-up, checked against real literature**: `ERBB3 ↔ ABL1` (CPTAC's strongest RTK-NRTK
correlation with no matching STRING edge) is a genuine, documented interaction, not a spurious
correlation. *The ERBB3 receptor in cancer and cancer gene therapy* (Cancer Gene Therapy,
Nature) states ERBB3 downstream signaling interacts with several partners including SRC, ABL,
rasGAP, and SYK — confirming real literature support for exactly the gap this cross-check was
built to surface: a real, biologically documented interaction that STRING's predicted network
simply doesn't capture.

## Real 90-kinase panel includes LMTK2/LMTK3, contradicting this project's own audit's explicit decision

A kinase-panel audit (`docs/kinase_panel_audit/kinase_audit_report.md`, dated March 2026, real
citations: Robinson 2000, Manning 2002, Trenker & Jura 2020/2021) corrected the original
54-RTK/29-NRTK (83 total) gene lists to 56 RTKs + 32 NRTKs = 88, adding 5 real, cited,
TNBC-relevant genes (STYK1, EPHA10, PTK6, SRMS, TXK) — and **explicitly excluded** LMTK2 and
LMTK3, stating plainly: "LMTK2 and LMTK3 are deliberately excluded from the RTK_GENES list in
this pipeline... This decision should be stated explicitly in the Methods section of any
resulting manuscript," since current classification (Trenker & Jura) reclassifies them as
serine/threonine kinases, not tyrosine kinases.

**Confirmed by an exact diff, not assumed**: `data/raw/kinases/kinase_90_list.txt` — the real
panel used for every CTS score, STRING centrality value, DGIdb query, and survival association
computed anywhere in this project — differs from the audit's validated 88-gene list by
*exactly* `LMTK2` and `LMTK3`, the two genes the audit says shouldn't be there. The audit's own
methodological decision was determined but never actually applied to the working panel.

**Not yet resolved.** Two of ninety kinases is a small fraction, and if their CTS/PairCTS
scores are unremarkable this may not meaningfully change any headline result — but that hasn't
been checked. Before publishing anything from this project, either (a) confirm LMTK2/LMTK3
don't materially affect any reported ranking and note their presence as a stated deviation from
the audit, or (b) regenerate `kinase_90_list.txt` and every downstream CTS/PairCTS/TripletCTS
result from the audited 88-gene list. Also note: a separate planning document
(`RTK_NRTK_Redundancy_Functions.md`) still cites the old, pre-audit 58-RTK count
(58×57/2 = 1,653 RTK-RTK pairs) — that reference material needs updating to 56 too, independent
of the pipeline-panel question above.

## HCOS pathway-overlap bias, self-identified in the authoritative grant proposal

`mc-ore-prototype/docs/Grant_Proposal.pdf` (v1.1) identifies a specific, real scoring bias in
the current HCOS formula: for the confirmed focal patient (TCGA-AO-A128), a second,
independently-generated hypothesis (`adavosertib + venetoclax`, via NMD-aware TP53 zygosity
resolution) is mechanistically sound but ranks lower than the headline
`afatinib+alpelisib+trastuzumab` triplet — not because it's biologically worse, but because
HCOS rewards shared-pathway overlap, and TP53's pathway doesn't overlap with the other three
genes'. The proposal itself proposes a pre-registered wet-lab test for this (Aim 3): if the
lower-HCOS regimen matches or beats the higher-HCOS one in real synergy assays, that's direct
evidence the pathway-overlap term needs down-weighting relative to evidence-strength and
mechanistic-complementarity terms. Not something to fix in code without that wet-lab evidence
first — noted here as a real, stated limitation, not a bug.

## Not clinically validated

Every result across every track in this repository is a computational research hypothesis
intended to prioritize candidates for wet-lab testing — not a clinical treatment
recommendation, and not validated by experimental synergy assay, clinical trial, or physician
oversight.
