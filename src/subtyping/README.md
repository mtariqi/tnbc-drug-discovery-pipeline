# subtyping/

## Built and verified: tnbc_de_novo_subtyping.py

Track C: TNBC molecular subtyping (Lehmann/Burstein-style: BL1, BL2, M, MSL, LAR, IM), via NMF
clustering on the confirmed-TNBC expression subset (reuses `restrict_to_tnbc()` from
`src/survival/tnbc_specific_survival.py`) plus marker-gene characterization. Necessarily
restricted to TNBC-only patients (~150-200), unlike Track A, since these subtypes aren't
defined outside TNBC.

**Real bug found and fixed during smoke-testing, not assumed correct:** the initial
argmax-over-raw-NMF-loadings cluster assignment silently failed for one of three synthetic
test groups — different NMF components can have very different natural scales (e.g. one
subtype's component ~1.2-1.5, another's genuinely elevated but only ~0.1-0.2), so raw argmax
unfairly favored whichever component happened to have the largest scale. Fixed by z-scoring
each component's loadings across patients before taking argmax. Confirmed: recovered 100%
cluster purity against known synthetic ground truth after the fix (was 40% purity for one
group before it).

**Honest limitation, stated in the module's own docstring:** there is no Lehmann/Burstein
ground-truth label in the TCGA clinical file to validate real cluster calls against — the TCGA
clinical file only has PAM50 (Basal/LumA/LumB/Her2/Normal), not TNBC-intrinsic subtypes. Marker
gene characterization is a defensible inference from expression pattern, not a confirmed
classification. Every real run should report the silhouette score and full marker-enrichment
table alongside cluster labels, not just the labels themselves.

## Next step

Run `run_tnbc_subtyping_pipeline()` against your real TNBC-restricted expression data
(the same `~/rtk_nrtk_tnbc/data/processed/tcga_brca/` files Track A uses, filtered via
`restrict_to_tnbc()`) to get real cluster assignments — this has only been verified against
synthetic data with known ground truth so far, same status Track A was in before its real-data
run.

