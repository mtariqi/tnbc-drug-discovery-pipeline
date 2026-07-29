"""
Resolve drug names -> SMILES strings, the missing join step flagged in
real_drugcomb_loader.py's docstring.

ORDER OF OPERATIONS, cheapest/most-reliable first:
  1. Check whether DrugComb's own bulk download already includes a drug metadata
     file with SMILES (many drug-combination portals ship one alongside the
     summary file, commonly named "drugs.csv" or similar). If you have this,
     USE IT INSTEAD of the PubChem fallback below -- it avoids a name-matching
     problem entirely, since it's already keyed to DrugComb's own drug IDs/names.
     Run inspect_drugcomb_csv() from real_drugcomb_loader.py on that file first to
     confirm its real column names (likely something like "dname"/"smiles" or
     "cid"/"smiles" -- unverified, check before trusting).
  2. For any drug NOT covered by DrugComb's own metadata, this module queries
     PubChem's public PUG REST API by name to fetch an Isomeric SMILES.

IMPORTANT CAVEATS, stated plainly:
  - This makes live network requests. It was written but NOT executed this session
    (no network access in the environment this was built in) -- run it yourself
    and treat the first real run as the actual test, same as other real-data
    loaders in this pipeline.
  - Biologics (e.g. trastuzumab, pembrolizumab) will typically return NO result
    from a small-molecule SMILES lookup, or an unreliable one for very large
    structures. This module does not guess -- a failed or suspiciously-large
    lookup is reported as a gap for you to resolve manually (drop the drug, or
    add a fallback ID-embedding path per synergy_gnn.py's docstring), not
    silently papered over.
  - PubChem name matching is not perfect: trade names, salt forms, and generic
    names can all fail to match. Check the "not found" list against your actual
    drug list by hand before assuming those drugs are unavailable in PubChem.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import requests

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
MAX_REASONABLE_ATOMS = 100  # a hand-picked ceiling; a lookup returning a SMILES this
                             # large for a "drug" is suspicious and worth a manual look
                             # rather than blind trust (could be a biologic's rare
                             # small-molecule-annotated fragment, not the real active drug).


def fetch_smiles_from_pubchem(drug_name: str, retry: int = 3, sleep_seconds: float = 0.3) -> str | None:
    """Single-drug PubChem lookup by name -> Isomeric SMILES. Returns None (not an
    exception) on any failure, so a caller processing a whole drug list doesn't
    stop at the first miss -- see resolve_drug_list() for how misses are reported.
    """
    url = f"{PUBCHEM_BASE}/compound/name/{drug_name}/property/IsomericSMILES/JSON"
    for attempt in range(retry):
        try:
            response = requests.get(url, timeout=10)
            if response.status_code == 200:
                data = response.json()
                return data["PropertyTable"]["Properties"][0]["IsomericSMILES"]
            if response.status_code == 404:
                return None  # genuinely not found, not a transient error -- don't retry
        except (requests.RequestException, KeyError, IndexError, ValueError):
            pass
        time.sleep(sleep_seconds * (attempt + 1))  # simple backoff
    return None


def resolve_drug_list(
    drug_names: list[str],
    cache_path: str | Path = "data/processed/drug_smiles_cache.json",
    known_smiles: dict[str, str] | None = None,
) -> dict[str, str]:
    """Resolve every drug name in `drug_names` to a SMILES string, checking
    `known_smiles` (e.g. parsed from DrugComb's own metadata file) first, then a
    local cache (so re-running doesn't re-hit PubChem for drugs already resolved),
    then falling back to a live PubChem query for anything still missing.

    Returns {drug_name: smiles} for every drug that resolved. Drugs that did not
    resolve are printed explicitly, not silently dropped -- decide what to do with
    them (exclude, or investigate the name mismatch) before proceeding.
    """
    known_smiles = known_smiles or {}
    cache_path = Path(cache_path)
    cache: dict[str, str] = {}
    if cache_path.exists():
        cache = json.loads(cache_path.read_text())
        print(f"Loaded {len(cache)} cached SMILES from {cache_path}.")

    resolved: dict[str, str] = {}
    unresolved: list[str] = []
    newly_fetched = 0

    for name in drug_names:
        if name in known_smiles:
            resolved[name] = known_smiles[name]
        elif name in cache:
            resolved[name] = cache[name]
        else:
            smiles = fetch_smiles_from_pubchem(name)
            if smiles is not None:
                resolved[name] = smiles
                cache[name] = smiles
                newly_fetched += 1
            else:
                unresolved.append(name)

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(json.dumps(cache, indent=2))

    print(f"Resolved {len(resolved)} of {len(drug_names)} drugs "
          f"({newly_fetched} newly fetched from PubChem this run, cache updated).")
    if unresolved:
        print(f"UNRESOLVED ({len(unresolved)}) -- check these by hand; likely biologics, "
              f"trade names, or salt-form mismatches, not necessarily unavailable in "
              f"PubChem under a different name: {unresolved}")
    return resolved


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--drug-list", required=True, help="Path to a text file, one drug name per line.")
    parser.add_argument("--cache", default="data/processed/drug_smiles_cache.json")
    parser.add_argument("--output", default="data/processed/drug_smiles_resolved.json")
    args = parser.parse_args()

    with open(args.drug_list) as f:
        drug_names = [line.strip() for line in f if line.strip()]

    resolved = resolve_drug_list(drug_names, cache_path=args.cache)
    Path(args.output).write_text(json.dumps(resolved, indent=2))
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
