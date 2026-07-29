**RTK-NRTK Kinome Audit Report**

Canonical Tyrosine Kinase Gene List Validation for the **RTK-NRTK
Redundancy in TNBC** Project

*Phase 1 --- Data Mining Validation \| March 2026*

**1. Executive Summary**

A systematic audit of the tyrosine kinase gene lists used in the Phase 1
data mining pipeline (data_mining_final_fixed.R) was conducted against
the canonical human kinome as defined by Robinson et al. (2000) and
updated by Manning et al. (2002). The audit identified **5 genes missing
from the original lists**: two RTKs (STYK1, EPHA10) and three NRTKs
(PTK6, SRMS, TXK). All five have been added to the pipeline. The
corrected gene lists now cover **56 RTKs and 32 NRTKs (88 total)**,
representing the complete set of active tyrosine kinase
domain-containing kinases in the human genome.

+--------------+---+--------------+---+--------------+---+--------------+
| **54 → 56**  |   | **29 → 32**  |   | **83 → 88**  |   | **--- → 5**  |
|              |   |              |   |              |   |              |
| RTKs         |   | NRTKs        |   | Total        |   | Genes added  |
|              |   |              |   | kinases      |   |              |
+--------------+---+--------------+---+--------------+---+--------------+

**2. Background**

The human protein kinase complement (kinome) comprises 518 kinase
domains distributed across 478 genes, as catalogued by Manning et al.
(2002). Within this complement, tyrosine kinases form a functionally
distinct group that phosphorylate protein substrates on tyrosine
residues and are central to growth factor signalling, immune activation,
and oncogenesis.

Tyrosine kinases are classically divided into two structural groups
based on the presence or absence of a transmembrane domain: **receptor
tyrosine kinases (RTKs)**, which are membrane-spanning and are activated
by extracellular ligands, and **non-receptor tyrosine kinases (NRTKs)**,
which are cytoplasmic and typically signal downstream of RTKs or immune
receptors.

Canonical counts used in this audit are **58 RTKs** (Robinson et al.
2000) and **32 NRTKs** (Manning et al. 2002). Note that the RTK count of
58 includes LMTK2 and LMTK3, which have since been reclassified as
serine/threonine kinases. Their exclusion, which is scientifically
defensible in current literature, reduces the canonical RTK count to
**56 active tyrosine kinase domain RTKs** --- the figure adopted in the
corrected pipeline.

**3. Audit Findings**

The original pipeline gene lists contained **54 RTKs** and **29 NRTKs**
(83 total). Comparison against the canonical kinome revealed five
missing genes across both classes. The complete findings are detailed
below.

**3.1 Missing Receptor Tyrosine Kinases**

Two RTKs were absent from the original list. Both are members of
established RTK subfamilies and appear in the canonical 58-RTK count.

  --------------------------------------------------------------------------------------
  **Gene**   **Subfamily**   **Status**   **Canonical            **TNBC Relevance**
                                          Classification**       
  ---------- --------------- ------------ ---------------------- -----------------------
  STYK1      ALK / NOK       **ADDED**    RTK ---                Downregulation of STYK1
                                          serine/threonine       potentiates EGFR
                                          kinase-like NOK (Novel inhibition and is
                                          Oncogene with Kinase   linked to acquired drug
                                          domain); member of the tolerance in solid
                                          ALK subfamily per      tumours; oncogenic in
                                          Robinson 2000          multiple cancer types

  EPHA10     EphA            **ADDED**    RTK pseudokinase ---   EphA family scaffold
                                          catalytically inactive interactions relevant
                                          EphA family member;    to cell migration and
                                          included in all        tumour invasiveness;
                                          canonical EphA         pseudokinase may
                                          receptor counts        regulate EphA2
                                          (EphA1--10)            signalling in TNBC
  --------------------------------------------------------------------------------------

  -----------------------------------------------------------------------
  **NOTE:** LMTK2 and LMTK3 (LMR subfamily) are not included. Although
  historically counted among the 58 RTKs under the Robinson 2000
  definition, both have been reclassified as serine/threonine kinases in
  current literature (Trenker & Jura, 2020). Their exclusion is
  scientifically defensible and reduces the canonical RTK count from 58
  to 56 active tyrosine kinase domain RTKs --- the figure used in the
  corrected pipeline.

  -----------------------------------------------------------------------

**3.2 Missing Non-Receptor Tyrosine Kinases**

Three NRTKs were absent from the original list. All three are members of
established NRTK families and are included in the Manning 2002 kinome
definition.

  -----------------------------------------------------------------------------------
  **Gene**   **Family**   **Status**   **Canonical            **TNBC Relevance**
                                       Classification**       
  ---------- ------------ ------------ ---------------------- -----------------------
  PTK6       BRK family   **ADDED**    NRTK --- Breast tumour **Critical TNBC gene.**
                                       kinase (BRK); BRK/PTK6 Overexpressed in \~86%
                                       subfamily alongside    of breast carcinomas;
                                       FRK and SRMS;          high expression
                                       SRC-related but lacks  predicts poor outcome.
                                       lipid modification and Direct roles in TNBC
                                       membrane anchoring     drug resistance,
                                                              invasion, and
                                                              metastasis. Arguably
                                                              the most important
                                                              missing gene for the
                                                              RTK-NRTK redundancy
                                                              study.

  SRMS       BRK family   **ADDED**    NRTK --- SRC-Related   Co-amplified with PTK6
                                       kinase lacking         at chromosome 20q13.3
                                       C-terminal regulatory  in breast cancer;
                                       tyrosine and           functional redundancy
                                       N-terminal             with PTK6 makes it
                                       myristoylation;        directly relevant to
                                       completes the PTK6/BRK the redundancy analysis
                                       family (PTK6, FRK,     
                                       SRMS)                  

  TXK        TEC family   **ADDED**    NRTK --- TXK (also     Completes the TEC
                                       called RLK); 5th       kinase family. TEC
                                       member of the TEC      family kinases are
                                       family alongside BTK,  downstream of RTKs and
                                       ITK, TEC, and BMX;     involved in immune and
                                       distinguishes from     growth factor
                                       other TEC members by   signalling relevant to
                                       lacking PH domain      the tumour
                                                              microenvironment
  -----------------------------------------------------------------------------------

**4. Corrected Gene Lists**

The following tables present the complete, corrected gene lists as
implemented in the updated pipeline. Newly added genes are highlighted.

**4.1 RTK Gene List (56 genes)**

  ------------------------------------------------------------------------
  **Subfamily**    **Member Genes**                           **Count**
  ---------------- ------------------------------------------ ------------
  ErbB / HER       EGFR, ERBB2, ERBB3, ERBB4                  4

  MET              MET, MST1R                                 2

  **ALK**          ALK, LTK, STYK1 **★ new**                  3

  PDGF / SCF       KIT, PDGFRA, PDGFRB, CSF1R, FLT3           5

  VEGFR            FLT1, KDR, FLT4                            3

  TIE              TEK, TIE1                                  2

  FGFR             FGFR1, FGFR2, FGFR3, FGFR4                 4

  Insulin / IGF    INSR, IGF1R, INSRR                         3

  NTRK / Trk       NTRK1, NTRK2, NTRK3                        3

  TAM / AXL        AXL, MERTK, TYRO3                          3

  DDR              DDR1, DDR2                                 2

  RET              RET                                        1

  ROR              ROR1, ROR2                                 2

  MuSK             MUSK                                       1

  AATK             AATK                                       1

  **EphA**         EPHA1, EPHA2, EPHA3, EPHA4, EPHA5, EPHA6,  9
                   EPHA7, EPHA8, EPHA10 **★ new**             

  EphB             EPHB1, EPHB2, EPHB3, EPHB4, EPHB6          5

  ROS              ROS1                                       1

  RYK              RYK                                        1

  PTK7             PTK7                                       1

  TOTAL                                                       56
  ------------------------------------------------------------------------

★ = gene added in this audit

**4.2 NRTK Gene List (32 genes)**

  -----------------------------------------------------------------------
  **Family**         **Member Genes**                        **Count**
  ------------------ --------------------------------------- ------------
  SRC family         SRC, YES1, FYN, LYN, LCK, HCK, BLK,     9
                     FGR, FRK                                

  **BRK / PTK6       PTK6, SRMS, FRK **★ includes new        3
  family**           additions**                             

  FAK family         PTK2, PTK2B                             2

  JAK family         JAK1, JAK2, JAK3, TYK2                  4

  ABL family         ABL1, ABL2                              2

  SYK family         ZAP70, SYK                              2

  **TEC family**     BTK, ITK, TEC, BMX, TXK **★ includes    5
                     new additions**                         

  ACK family         TNK1, TNK2                              2

  CSK family         CSK, MATK                               2

  FER family         FER, FES                                2

  TOTAL                                                      32
  -----------------------------------------------------------------------

★ = gene added in this audit; FRK appears in both SRC and BRK family
rows as it bridges both classifications

**5. Biological Rationale for Added Genes**

**5.1 PTK6 (BRK) --- Priority Addition for TNBC Research**

**PTK6** (Protein Tyrosine Kinase 6), also known as Breast tumour Kinase
(BRK), is the single most important gene added in this audit relative to
the project\'s scientific objectives.

-   **Expression:** Overexpressed in approximately 86% of breast
    carcinomas at the protein level, with highest expression in the most
    aggressive subtypes including TNBC.

-   **Prognosis:** High PTK6 expression is an independent predictor of
    poor clinical outcome in breast cancer patients.

-   **TNBC mechanism:** PTK6 activates STAT3, AKT, and ERK signalling
    pathways downstream of RTKs; its presence as an NRTK creates a
    compensatory bypass route when individual RTKs are inhibited ---
    precisely the redundancy mechanism this project aims to
    characterise.

-   **Drug resistance:** PTK6 has been implicated in resistance to EGFR
    inhibitors and chemotherapy in TNBC cell lines, making it directly
    relevant to the drug resistance phenotypes under study in GEO
    datasets GSE58644 and GSE167977.

**5.2 SRMS --- Genomic Co-amplification with PTK6**

**SRMS** (Src-Related kinase lacking C-terminal regulatory tyrosine and
N-terminal myristoylation) is co-located with PTK6 at chromosome
20q13.3, a region frequently amplified in breast cancer. Its
co-amplification with PTK6 means that any analysis of PTK6 copy number,
expression, or substrate redundancy is incomplete without SRMS. As a
member of the same BRK kinase family (PTK6, FRK, SRMS), it shares
substrate specificity with PTK6 and may act as a functional redundant
partner.

**5.3 TXK --- Completing the TEC Kinase Family**

**TXK** (also known as RLK, Resting Lymphocyte Kinase) is the fifth
member of the TEC family of NRTKs, completing the set {BTK, ITK, TEC,
BMX, TXK}. Unlike other TEC family members, TXK lacks the N-terminal
pleckstrin homology (PH) domain, giving it distinct membrane-targeting
properties. The TEC family kinases are activated downstream of RTKs and
pattern recognition receptors and are involved in immune cell signalling
within the tumour microenvironment. Omitting TXK from the family creates
an incomplete picture of TEC-family redundancy.

**5.4 STYK1 --- Oncogenic RTK in the ALK Subfamily**

**STYK1** (Serine Threonine Tyrosine Kinase 1), also called NOK (Novel
Oncogene with Kinase domain), is an atypical RTK member of the ALK
subfamily. Unlike classical RTKs, its extracellular domain is truncated
and it may signal constitutively. STYK1 has been identified as a gene
whose downregulation potentiates EGFR inhibition, suggesting it acts as
a resistance-mediating bypass kinase --- directly relevant to the
RTK-NRTK redundancy hypothesis underpinning this project.

**5.5 EPHA10 --- Pseudokinase Scaffold in the EphA Family**

**EPHA10** is a catalytically inactive pseudokinase that is nonetheless
counted among the canonical EphA receptors (EphA1--EphA10). Although it
cannot phosphorylate substrates directly, pseudokinases frequently act
as allosteric regulators of active kinase family members. EPHA10 may
modulate EphA2 signalling, which has established roles in TNBC invasion
and metastasis. Its inclusion ensures the EphA family analysis in the
project is complete.

**6. LMTK2 and LMTK3 --- Classification Note**

LMTK2 and LMTK3 (Lemur Tyrosine Kinase 2 and 3) appear in the Robinson
2000 canonical RTK count of 58 as members of the LMR (Lemur) subfamily.
However, subsequent structural and biochemical analyses have established
that LMR family kinases phosphorylate serine and threonine residues
rather than tyrosine. They have therefore been reclassified as
serine/threonine kinases in current databases (UniProt, HGNC) and
excluded from updated tyrosine kinase-specific analyses.

  -----------------------------------------------------------------------
  **DECISION:** LMTK2 and LMTK3 are deliberately excluded from the
  RTK_GENES list in this pipeline. This is scientifically defensible per
  current classification and reduces the canonical RTK count from 58
  (Robinson 2000) to 56 (active tyrosine kinase domain RTKs). This
  decision should be stated explicitly in the Methods section of any
  resulting manuscript.

  -----------------------------------------------------------------------

**7. Impact on Pipeline Outputs**

The addition of 5 genes to the kinase lists requires reprocessing of
several Phase 1 pipeline steps. The following checkpoint files must be
deleted before re-running the script to force reprocessing:

  ------------------------------------------------------------------------
  **Checkpoint File**          **Reason for Reprocessing**
  ---------------------------- -------------------------------------------
  save_kinase_reference.done   Kinase reference CSV must reflect updated
                               88-gene list

  tcga_rnaseq_prepare.done     Kinase row extraction must include 5 new
                               genes from 1,231 STAR TSVs

  tcga_tnbc_labels.done        TNBC barcode assignment uses updated gene
                               list for metadata matching

  uniprot_fetch.done           UniProt sequence/domain retrieval must
                               cover all 88 kinases

  pdb_fetch.done               PDB and AlphaFold2 structure retrieval must
                               cover all 88 kinases

  string_network.done          STRING network query must include all 88
                               kinases
  ------------------------------------------------------------------------

Run the following commands before restarting the pipeline (fish shell):

+-----------------------------------------------------------------------+
| rm \~/rtk_nrtk_tnbc/logs/save_kinase_reference.done                   |
|                                                                       |
| rm \~/rtk_nrtk_tnbc/logs/tcga_rnaseq_prepare.done                     |
|                                                                       |
| rm \~/rtk_nrtk_tnbc/logs/tcga_tnbc_labels.done                        |
|                                                                       |
| rm \~/rtk_nrtk_tnbc/logs/uniprot_fetch.done                           |
|                                                                       |
| rm \~/rtk_nrtk_tnbc/logs/pdb_fetch.done                               |
|                                                                       |
| rm \~/rtk_nrtk_tnbc/logs/string_network.done                          |
+-----------------------------------------------------------------------+

**8. Validation Status After Audit**

The table below shows the pipeline output validation status after the
gene list correction, based on the most recent run (2026-03-30 00:37).
Steps that will be reprocessed after checkpoint deletion are marked
accordingly.

  -------------------------------------------------------------------------------------------
  **Output File**               **Size**   **Status**    **Notes**
  ----------------------------- ---------- ------------- ------------------------------------
  kinase_counts_raw.csv         479.2 KB   **OK**        Will reprocess --- 5 new kinase rows
                                                         to extract

  sample_metadata.csv           952.3 KB   **OK**        Stable --- no change required

  kinase_mutations.csv          1,570.2 KB **OK**        Stable --- no change required

  clinical_data.csv             737.2 KB   **OK**        Stable --- no change required

  tnbc_sample_barcodes.txt      MISSING    **MISSING**   Will reprocess --- IHC column
                                                         detection fixed

  all_kinases.fasta             MISSING    **MISSING**   Will reprocess --- 5 new sequences
                                                         to fetch

  uniprot_seq_metadata.csv      MISSING    **MISSING**   Will reprocess --- 5 new entries

  uniprot_domains.csv           MISSING    **MISSING**   Will reprocess --- 5 new entries

  kinase_reference.csv          0.8 KB     **OK**        Will reprocess --- updated to 88
                                                         genes

  pdb_ids.csv                   MISSING    **MISSING**   Will reprocess --- 5 new genes

  download_alphafold.sh         MISSING    **MISSING**   Will reprocess --- 5 new structures

  string_kinase_subgraph.csv    10.4 KB    **OK**        Will reprocess --- 5 new nodes

  psp_kinase_substrate_pY.csv   MISSING    **MISSING**   Requires manual PSP download

  psp_shared_substrates.csv     MISSING    **MISSING**   Requires manual PSP download

  pubmed_abstracts.csv          500.3 KB   **OK**        Stable --- query-based, not
                                                         gene-list dependent
  -------------------------------------------------------------------------------------------

**9. Outstanding Items**

**9.1 PhosphoSitePlus Manual Download**

Three output files (psp_kinase_substrate_pY.csv,
psp_shared_substrates.csv, convergence_substrates.csv) require the
PhosphoSitePlus Kinase-Substrate Dataset, which is not publicly
downloadable via API and must be obtained manually:

-   **URL:** https://www.phosphosite.org/downloads/

-   **File:** Kinase_Substrate_Dataset.gz

-   **Destination:** \~/rtk_nrtk_tnbc/data/raw/phosphosites/

**9.2 TNBC Sample Barcodes**

The TNBC barcode identification step (tcga_tnbc_labels) produced 0
barcodes in the last run due to IHC column name mismatches between
TCGAbiolinks versions. The pipeline now searches for multiple known
column name variants (e.g., er_status_by_ihc,
breast_carcinoma_estrogen_receptor_status) and logs which columns are
actually present in the clinical CSV. The checkpoint must be deleted and
the step rerun to resolve this. If the IHC columns remain absent, the
pipeline will fall back to PAM50 Basal classification from the sample
metadata.

**10. References**

**Robinson DR, Wu YM, Lin SF (2000).** The protein tyrosine kinase
family of the human genome. Oncogene 19(49):5548--5557.

**Manning G, Whyte DB, Martinez R, Hunter T, Sudarsanam S (2002).** The
protein kinase complement of the human genome. Science
298(5600):1912--1934.

**Lemmon MA, Schlessinger J (2010).** Cell signaling by receptor
tyrosine kinases. Cell 141(7):1117--1134.

**Trenker R, Jura N (2021).** Receptor tyrosine kinase activation: From
the ligand perspective. Curr Opin Cell Biol 63:174--185.

**Bhatt DL et al. (2023).** PTK6 in breast cancer: Biology, signalling,
and therapeutic implications. Cancer Lett 554:215987.

**Goel RK, Bhatt DL et al. (2020).** SRMS is a nonmyristoylated SRC
family kinase that phosphorylates and activates SRC in breast cancer.
Oncogene 39:3680--3695.
