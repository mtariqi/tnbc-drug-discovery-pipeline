import sys
sys.path.insert(0, '.')

from mdcoe import GENE_DRUGS, resolve_tp53_drugs
from build_real_gene_drugs import build_real_gene_drugs
from agentic_regimen_discovery import run_discovery_loop, merge_into_candidate_pool

genes = ["EGFR", "PTEN", "FLT1", "TP53", "ERBB2", "TYK2"]
kinase_panel = [k.strip() for k in open("/home/mtariq/rtk_nrtk_tnbc/data/raw/kinases/kinase_90_list.txt")]
dgidb_path = "/home/mtariq/rtk_nrtk_tnbc/data/processed/dgidb/dgidb_interactions.tsv"

# Build the correct curated baseline -- TP53 needs resolve_tp53_drugs(),
# not the raw GENE_DRUGS entry, which is generic and alteration-agnostic
curated_gene_drugs = dict(GENE_DRUGS)
curated_gene_drugs["TP53"] = resolve_tp53_drugs("Nonsense_Mutation")  # confirm this alteration type against the real MAF file

real_gene_drugs = build_real_gene_drugs(dgidb_path, genes, all_kinase_panel=kinase_panel)
discovery = run_discovery_loop(genes, curated_gene_drugs, real_gene_drugs, min_candidates=4)

if not discovery.empty:
    print(discovery[["gene", "candidate_drug", "pmid", "dgidb_confirmed", "mention_context"]])
    expanded_pool = merge_into_candidate_pool(discovery, curated_gene_drugs)
    print(expanded_pool["TP53"])
else:
    print("No candidates found.")
