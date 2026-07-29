# scoring/

## Real code (not placeholders)

- `kinase_scoring_pipeline.py` — `DrugScore`, `DrugPairScore`, `TripletCTS` and friends.
  **Verified in this session**: its own synthetic demo runs correctly (ranks afatinib highest
  among the 6 demo drugs, produces sensible pair/triplet scores) — but this is still synthetic
  demo data, not yet run against the real `~/rtk_nrtk_tnbc/data/` files.
- `validate_pipeline_inputs.py` — the 11-point pre-flight check. **Verified against the real
  data**: 16 passed, 2 warnings (STRING network coverage 85/90, FAERS toxicity coverage
  95/117 — see limitations.md), 0 failures.
- `run_pairs_and_triplets.py` — runs real PairCTS/TripletCTS against the real, already-computed
  `CTS_all_90_kinases.tsv`. **Run for real against the actual data this session, with three
  real bugs found and fixed along the way** (not assumed correct from the code alone):
  1. DGIdb's raw uppercase/salt-form drug names (`ERLOTINIB HYDROCHLORIDE`) didn't match the
     lowercase plain-generic names in `drug_list.txt`/`GENE_DRUGS` — 0/117 drugs matched
     before the fix.
  2. Generic mechanism/class labels (`antineoplastic agent`, `btk inhibitor`) were present in
     DGIdb's raw output as if they were specific drug names — an exact-phrase stoplist now
     rejects these, deliberately NOT a token-level rule (which would have also incorrectly
     rejected real compound names like `paclitaxel protein-bound`).
  3. `drug_list.txt` itself has 3 inconsistent entries carrying their salt-form suffix as part
     of the curated name (`cabozantinib s-malate`, `dacomitinib anhydrous`,
     `dasatinib anhydrous`) — fixed by normalizing both sides of the match identically.
  Final real result: **116/117 drugs scoreable** — the sole holdout (`antineoplastic agent`)
  is a genuine curation issue in `drug_list.txt` itself (worth removing from that file
  directly), not a code problem.

## Confirmed real finding, not yet resolved

Running this for real surfaced a genuine methodological limitation, not a bug: `PairCTS`'s
`best_target_pair_cts()` takes the max across all target-pair combinations, which produces
heavy score ties (`0.656483` for most top-20 real pairs) whenever multiple drug pairs collapse
onto the same dominant kinase pair (confirmed: two independent drug pairs both resolved to the
literal same winning kinase pair, `ERBB2`+`EPHA5`). See `docs/limitations.md` for the full
writeup and candidate fixes — this is a real scoring-design decision worth making deliberately,
not silently patching.

## PairCTS ceiling-effect: resolved (partially), adopted as default

Two remediation attempts, both tested against the real 116-drug data (see `docs/limitations.md`
for full detail): the first (`primary_target_pair_scoring.py`) was rejected — it made things
worse. The second (`redundancy_penalized_pair_scoring.py`) gave a genuine, if partial,
improvement for top-N reporting specifically (top-20 ties: 20/20 → 7/20; distinct winning
kinase pairs in top 20: ~1 → 13), at the cost of the full-ranking's aggregate uniqueness. **This
has been adopted as `run_pairs_and_triplets.py`'s default** for its printed top-N pairs and
triplet-building drug pool (weight=0.15) — the original, unpenalized `rank_all_pairs()` remains
available and unchanged in `kinase_scoring_pipeline.py` for anyone who wants the full ranking.

## Next step

Run `kinase_scoring_pipeline.py`'s `compute_cts()` against the real (non-synthetic)
`~/rtk_nrtk_tnbc/data/` inputs directly — `run_pairs_and_triplets.py` above loads an
*already-computed* `CTS_all_90_kinases.tsv`, so the CTS-computation step itself is still only
verified against synthetic demo data, not real data end-to-end.

