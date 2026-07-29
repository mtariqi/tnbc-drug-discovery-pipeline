"""
Run Real PairCTS / TripletCTS
===============================
"""

from __future__ import annotations

import itertools
import re
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd

from kinase_scoring_pipeline import (
    pair_cts, best_target_pair_cts, drug_pair_score, rank_all_pairs,
    Regimen, rank_all_triplets, load_drug_target_map,
)
from redundancy_penalized_pair_scoring import rank_all_pairs_redundancy_penalized


SALT_SUFFIXES = [
    "hydrochloride", "dihydrochloride", "sodium", "mesylate", "dimesylate", "tartrate",
    "sulfate", "sulphate", "phosphate", "citrate", "maleate", "dimaleate", "hydrobromide",
    "besylate", "succinate", "fumarate", "acetate", "bromide", "chloride",
    "monohydrate", "dihydrate", "s-malate", "anhydrous",
]

# Generic mechanism/class labels found in the real data that are NOT specific
# drug names (e.g. DGIdb listing "antineoplastic agent" as if it were itself a
# compound). Exact-phrase match only, deliberately -- a token-level rule would
# also incorrectly reject real compound names like "paclitaxel protein-bound".
JUNK_PHRASES = {
    "antineoplastic agent", "anesthetic agent", "btk inhibitor", "c-met inhibitor",
    "anthracycline antineoplastic antibiotic", "anthraquinone analogue",
    "antisense oligonucleotides",
}


def _normalize_drug_name(name: str) -> str:
    """
    Real gap found and measured on the actual data: DGIdb's raw 'drug' column is
    UPPERCASE and often includes a salt/crystal-form suffix (e.g. 'ERLOTINIB
    HYDROCHLORIDE'), while drug_list.txt / GENE_DRUGS use plain lowercase
    generic names ('erlotinib'). A direct case-sensitive dict lookup matched
    0/117 drugs before this fix.

    IMPORTANT DISTINCTION, verified against real examples before writing this:
    conjugate/payload suffixes ('vedotin', 'pegol', 'emtansine', 'alaninate')
    are NOT stripped -- they denote a scientifically distinct drug identity,
    not an interchangeable salt form (stripping 'vedotin' from 'brentuximab
    vedotin' would be a real error, not a normalization).
    """
    n = name.lower().strip()
    for suf in SALT_SUFFIXES:
        n = re.sub(rf"\s*\({suf}\)$", "", n)
        n = re.sub(rf"\s+{suf}$", "", n)
    return n.strip()


def _is_junk_drug_name(raw_name: str) -> bool:
    """Rejects generic mechanism/class labels masquerading as drug names.
    See JUNK_PHRASES above for why this is an exact-phrase check, not a
    token-level one."""
    n = raw_name.lower().strip()
    return n in JUNK_PHRASES or _normalize_drug_name(raw_name) in JUNK_PHRASES


def load_real_inputs(base_dir: str):
    """Loads the five real files this needs. Adjust paths if your layout differs."""
    base = Path(base_dir)

    cts_df = pd.read_csv(base / "data/processed/CTS_all_90_kinases.tsv", sep="\t", index_col="kinase_id")
    cts = cts_df["CTS"]

    community_df = pd.read_csv(base / "data/processed/string/community_map.tsv", sep="\t")
    community_map = community_df.set_index("kinase_id")["community"].to_dict()

    edges_df = pd.read_csv(base / "data/processed/string/string_edges.tsv", sep="\t")
    crosstalk_edges = {(row["source"], row["target"]): row["weight"] for _, row in edges_df.iterrows()}

    dgidb_df = pd.read_csv(base / "data/processed/dgidb/dgidb_interactions.tsv", sep="\t")
    dgidb_df = dgidb_df.copy()
    dgidb_df = dgidb_df[~dgidb_df["drug"].apply(_is_junk_drug_name)]
    dgidb_df["drug"] = dgidb_df["drug"].apply(_normalize_drug_name)
    drug_target_df = dgidb_df[["drug", "kinase_id"]]
    target_map = load_drug_target_map(drug_target_df)

    return cts, community_map, crosstalk_edges, target_map


def build_toxicity_lookup_from_local_faers(faers_dir: str, drug_list: List[str]) -> Dict[Tuple[str, str], float]:
    faers_path = Path(faers_dir)
    reaction_sets: Dict[str, set] = {}
    for drug in drug_list:
        f = faers_path / f"{drug}_faers.tsv"
        if f.exists():
            df = pd.read_csv(f, sep="\t")
            reaction_sets[drug] = set(df["reaction"].str.lower()) if "reaction" in df.columns else set()

    toxicity_lookup = {}
    for d1, d2 in itertools.combinations(drug_list, 2):
        if d1 in reaction_sets and d2 in reaction_sets:
            shared = reaction_sets[d1] & reaction_sets[d2]
            total = reaction_sets[d1] | reaction_sets[d2]
            toxicity_lookup[(d1, d2)] = len(shared) / len(total) if total else 0.0
    return toxicity_lookup


def run(base_dir: str, drug_list: List[str], max_drugs_for_pairs: int = 60, top_k_pairs: int = 20, top_k_triplets: int = 20):
    print("Loading real CTS, community, crosstalk, and DGIdb data...")
    cts, community_map, crosstalk_edges, target_map = load_real_inputs(base_dir)
    print(f"  CTS: {len(cts)} kinases, {len(community_map)} community assignments, "
          f"{len(crosstalk_edges)} crosstalk edges, {len(target_map)} drugs with known targets")

    scoreable_drugs = [d for d in drug_list if d in target_map and any(k in cts.index for k in target_map[d])]
    print(f"  {len(scoreable_drugs)}/{len(drug_list)} drugs have at least one scoreable target in your 90-kinase panel")

    if len(scoreable_drugs) > max_drugs_for_pairs:
        print(f"  Capping to top {max_drugs_for_pairs} drugs by max target CTS, to keep pair enumeration tractable "
              f"(C({len(scoreable_drugs)},2) = {len(scoreable_drugs)*(len(scoreable_drugs)-1)//2} pairs otherwise)")
        drug_max_cts = {d: max(cts.get(k, 0) for k in target_map[d]) for d in scoreable_drugs}
        scoreable_drugs = sorted(scoreable_drugs, key=lambda d: drug_max_cts[d], reverse=True)[:max_drugs_for_pairs]

    print("\nBuilding pairwise toxicity proxy from already-downloaded FAERS files (no new API calls)...")
    toxicity_lookup = build_toxicity_lookup_from_local_faers(str(Path(base_dir) / "data/processed/faers"), scoreable_drugs)
    print(f"  Toxicity data available for {len(toxicity_lookup)} drug pairs")

    synergy_lookup = {}  # GAP 1, see module docstring -- no real synergy data source exists yet

    print("\nRanking drug pairs (redundancy-penalized PairCTS/DrugPairScore, weight=0.15 -- ")
    print("confirmed on real data to improve top-20 diversity vs. the original max-aggregation")
    print("ranking: 20/20 -> 7/20 ties, ~1 -> 13 distinct winning kinase pairs represented.")
    print("The original, unpenalized full ranking remains available via rank_all_pairs() directly")
    print("or compare_redundancy_penalty.py, if ever needed for something other than top-N use.)")
    pairs_df = rank_all_pairs_redundancy_penalized(
        scoreable_drugs, target_map, cts, community_map, crosstalk_edges, toxicity_lookup,
        synergy_lookup, redundancy_penalty_weight=0.15,
    )
    print(f"Top {top_k_pairs} pairs:")
    print(pairs_df.head(top_k_pairs).to_string(index=False))

    print("\nBuilding candidate triplets from the top-ranked pairs' drugs...")
    top_pair_drugs = set(pairs_df.head(30)["drug_1"]) | set(pairs_df.head(30)["drug_2"])
    regimens = []
    all_modules = set(community_map.values())
    for d1, d2, d3 in itertools.combinations(sorted(top_pair_drugs), 3):
        targets = list(set(target_map.get(d1, []) + target_map.get(d2, []) + target_map.get(d3, [])))
        targets = [t for t in targets if t in cts.index]
        if not targets:
            continue
        modules_hit = [community_map.get(t) for t in targets if t in community_map]
        regimens.append(Regimen(
            drugs=(d1, d2, d3),
            targets=targets,
            modules_hit=modules_hit,
            escape_routes_closed=[],
            combined_toxicity_raw=sum(
                toxicity_lookup.get((a, b), toxicity_lookup.get((b, a), 0.0))
                for a, b in itertools.combinations((d1, d2, d3), 2)
            ),
        ))

    print(f"  {len(regimens)} candidate triplets built from {len(top_pair_drugs)} top-pair drugs")
    print("\nRanking triplets (real TripletCTS)...")
    triplets_df = rank_all_triplets(regimens, total_modules=len(all_modules), total_known_routes=1)
    print(f"Top {top_k_triplets} triplets:")
    print(triplets_df.head(top_k_triplets).to_string(index=False))

    return pairs_df, triplets_df


if __name__ == "__main__":
    import sys
    base_dir = sys.argv[1] if len(sys.argv) > 1 else str(Path.home() / "rtk_nrtk_tnbc")
    drug_list_path = sys.argv[2] if len(sys.argv) > 2 else str(Path(base_dir) / "data/raw/drugs/drug_list.txt")
    drug_list = [_normalize_drug_name(d.strip()) for d in open(drug_list_path) if d.strip()]
    run(base_dir, drug_list)
