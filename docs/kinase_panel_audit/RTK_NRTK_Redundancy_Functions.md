**KINASE REDUNDANCY FUNCTIONS**

**RTK--RTK \| NRTK--NRTK \| RTK--NRTK**

**in Triple-Negative Breast Cancer**

*Comprehensive Redundancy Function Reference*

Mechanisms · Computational Functions · Expected Outputs · Innovation
Pathway

**Overview: The Three-Head Redundancy Framework**

Kinase redundancy in TNBC operates at three distinct but interconnected
levels. Each level has its own molecular logic, computational
measurement strategy, and clinical implication. This document defines,
for each level, the precise redundancy functions implemented in the
pipeline --- what they compute, how they compute it, what they output,
and why they represent a novel scientific contribution.

  ---------------------------------------------------------------------------------
  **Head**   **Pair     **Core Redundancy Logic**        **No. of  **Clinical
             Type**                                      Pairs**   Relevance**
  ---------- ---------- -------------------------------- --------- ----------------
  HEAD 1     RTK -- RTK Two RTKs compensate when they    58×57/2 = When one RTK is
                        share: kinase domain identity,   1,653     inhibited (e.g.
                        overlapping ligand-activated     pairs     EGFR), a related
                        dimerisation,                              RTK (e.g. ERBB3)
                        trans-phosphorylation capacity,            transactivates
                        and convergent downstream                  the same
                        substrate repertoire                       downstream
                                                                   cascade;
                                                                   co-inhibition
                                                                   required

  HEAD 2     NRTK --    Two NRTKs compensate when they   32×31/2 = When SRC is
             NRTK       share: SH2 binding specificity   496 pairs inhibited, YES1
                        (same pY motif recognition), SH3           or LYN can
                        proline-rich ligand selectivity,           substitute at
                        scaffold co-occupancy, and                 identical pY
                        mutual exclusivity in                      substrates;
                        phosphoproteome activation                 functional
                                                                   redundancy is
                                                                   invisible to
                                                                   sequence-based
                                                                   drug design

  HEAD 3     RTK --     An RTK and NRTK compensate when  58×32 =   The most
             NRTK       they converge on shared          1,856     dangerous
                        phosphosubstrates via            pairs     cross-class
                        independent upstream routes, or            redundancy: EGFR
                        when the NRTK transactivates the           inhibition
                        RTK (or vice versa), sustaining            activates SRC
                        downstream flux through                    which
                        alternative entry points                   independently
                                                                   sustains AKT and
                                                                   ERK; requires
                                                                   cross-class
                                                                   co-inhibition
  ---------------------------------------------------------------------------------

+-----------+----------------------------------------------------------+
| **1**     | **RTK -- RTK Redundancy**                                |
|           |                                                          |
| HEAD      | *Intra-class compensation between receptor tyrosine      |
|           | kinases*                                                 |
+-----------+----------------------------------------------------------+

**Biological Rationale**

Receptor tyrosine kinases are single-pass transmembrane proteins that
activate upon extracellular ligand binding, triggering kinase domain
autophosphorylation and downstream signalling. RTK families (EGFR/ERBB,
FGFR, VEGFR, PDGFR, MET/RON, EPH, TAM, AXL) share conserved kinase
domain architecture. Within and across families, RTK redundancy arises
from four converging mechanisms:

-   Trans-phosphorylation: ErbB family members (EGFR, ERBB2, ERBB3,
    ERBB4) form obligate heterodimers; inhibiting one partner is rescued
    by an alternate dimer configuration that re-activates the same
    signalling node

-   Ligand promiscuity: Multiple RTKs bind shared or structurally
    similar ligands (e.g. HGF activates both MET and RON; FGF ligands
    activate all four FGFRs), so loss of one receptor does not deplete
    the ligand signal

-   Substrate convergence: Different RTKs phosphorylate the same adapter
    proteins (GAB1, IRS1, SHC1, GRB2) and therefore activate identical
    RAS-MAPK and PI3K-AKT cascades downstream

-   Compensatory transcriptional upregulation: Inhibition of one RTK
    triggers feedback-driven transcriptional upregulation of a family
    member --- documented for EGFR inhibition driving MET and FGFR1
    upregulation in TNBC

**RTK-RTK Redundancy Functions**

The following seven computational functions are implemented in the
pipeline to quantify RTK-RTK redundancy. Together they produce the
RTK-RTK Redundancy Score (RRS₁) for all 1,653 RTK pairs.

  -----------------------------------------------------------------------
  **F1.1 Kinase Domain Sequence Identity Score (KDSI)**

  **MECHANISM:** Pairwise percentage identity of kinase domain sequences
  extracted by HMMER (Pfam PF00069). Computed on MUSCLE multiple sequence
  alignment of all 58 RTK kinase domains. Identity calculated over the
  aligned region only (gapped positions excluded). High KDSI means both
  RTKs have near-identical catalytic machinery and are likely to
  phosphorylate the same substrate tyrosine motifs.

  **INPUTS:** 58 RTK canonical sequences from UniProt; HMMER-extracted
  kinase domains (PF00069); MUSCLE MSA

  **OUTPUTS:** 58×58 pairwise KDSI matrix (0--100%); pairs with KDSI
  \>40% flagged as \'high molecular redundancy potential\'; annotated
  phylogenetic tree showing KDSI-based clustering

  **INNOVATION:** First systematic KDSI matrix for all 58 human RTKs
  simultaneously. Prior studies computed pairwise identity for specific
  family pairs. The full 58×58 matrix reveals cross-family redundancy
  pairs not predicted by phylogeny alone.
  -----------------------------------------------------------------------

  -----------------------------------------------------------------------
  **F1.2 Gatekeeper Residue Identity Function (GRIF)**

  **MECHANISM:** Extracts the gatekeeper residue (the amino acid at the
  position equivalent to EGFR T790) from the kinase domain MSA for all 58
  RTKs. The gatekeeper controls access to the ATP-binding back pocket and
  is the primary determinant of Type II inhibitor selectivity. Two RTKs
  sharing a gatekeeper identity are cross-inhibited by the same clinical
  compounds, meaning inhibiting one may not fully block the other if
  cross-reactivity is insufficient.

  **INPUTS:** MUSCLE MSA output from F1.1; crystal structure coordinates
  for gatekeeper position validation (PDB/AlphaFold2)

  **OUTPUTS:** Gatekeeper residue table for all 58 RTKs; gatekeeper
  identity matrix (binary, 58×58); list of RTK pairs sharing identical
  gatekeeper --- these are inhibitor cross-reactivity candidates

  **INNOVATION:** Gatekeeper comparison across all 58 RTKs simultaneously
  has not been published. This function identifies RTK pairs that
  clinical inhibitors fail to distinguish --- directly predicting which
  resistance events are pharmacologically addressable.
  -----------------------------------------------------------------------

  -----------------------------------------------------------------------
  **F1.3 Structural Binding Pocket Similarity Function (SBPS)**

  **MECHANISM:** 3D structural superposition of kinase domains using
  ProDy (Python). Pairwise RMSD computed for all 58×58 RTK pairs after
  aligning to the activation loop reference frame. fpocket then analyses
  the ATP-binding pocket of each structure: volume, druggability score,
  hydrophobicity, and electrostatic surface. Tanimoto coefficient
  computed between fpocket descriptors for each pair. Low RMSD + high
  Tanimoto = structurally redundant binding pockets that will respond to
  the same inhibitor scaffold.

  **INPUTS:** PDB crystal structures + AlphaFold2 models for all 58 RTKs;
  ProDy superposition; fpocket analysis

  **OUTPUTS:** 58×58 RMSD matrix (Å); 58×58 pocket Tanimoto matrix; list
  of pairs with RMSD \<2.0Å AND Tanimoto \>0.70 (high structural
  redundancy); Figure 2b heatmap data

  **INNOVATION:** First simultaneous structural comparison of all 58
  human RTK binding pockets. Identifies structurally redundant pairs that
  would be co-inhibited by the same clinical kinase inhibitor --- a key
  factor in predicting which drug combinations are truly orthogonal.
  -----------------------------------------------------------------------

  -----------------------------------------------------------------------
  **F1.4 Shared Phosphosubstrate Overlap Function (SPOF-RTK)**

  **MECHANISM:** For each RTK pair (i,j), computes the Jaccard index of
  their phosphosubstrate sets using PhosphoSitePlus (PSP)
  kinase-substrate data. Substrate set S(i) = all proteins phosphorylated
  by RTK i at any pY site. Jaccard(i,j) = \|S(i) ∩ S(j)\| / \|S(i) ∪
  S(j)\|. Additionally computes the raw shared substrate count and
  identifies \'convergence substrates\' --- substrates phosphorylated by
  3+ RTKs simultaneously (these are the most clinically critical nodes).

  **INPUTS:** PhosphoSitePlus human kinase-substrate-pY relationships; 58
  RTK substrate sets extracted

  **OUTPUTS:** Jaccard matrix for all 1,653 RTK pairs; shared substrate
  count per pair; convergence substrate list (pY sites targeted by ≥3
  RTKs); bipartite RTK-substrate network for Phase 4

  **INNOVATION:** The convergence substrate concept is novel: substrates
  phosphorylated redundantly by 3+ RTKs cannot be blocked by inhibiting
  any single RTK. These substrates --- likely GAB1, IRS1, SHC1, GRB2 ---
  are proposed as the first rationally defined class of phosphoproteomic
  TNBC resistance biomarkers.
  -----------------------------------------------------------------------

  -----------------------------------------------------------------------
  **F1.5 RTK Co-Expression Rewiring Function (CERF-RTK)**

  **MECHANISM:** Computes Spearman correlation of expression across
  TCGA-BRCA TNBC samples (VST-normalised) for all 1,653 RTK pairs.
  Additionally implements a \'compensatory shift score\': using GEO
  drug-resistance datasets, measures whether the expression of RTK j
  increases when RTK i is pharmacologically suppressed (treated vs
  baseline delta). High positive correlation + positive compensatory
  shift = strongest functional evidence that j compensates for i in real
  tumour tissue.

  **INPUTS:** TCGA-BRCA VST counts (TNBC samples only); GEO resistance
  datasets (GSE58644, GSE86948, GSE167977); DESeq2 DE results from Phase
  3

  **OUTPUTS:** Spearman matrix for 1,653 RTK pairs across TNBC;
  compensatory shift scores per pair; top 20 RTK pairs with highest
  combined co-expression + compensatory shift evidence; WGCNA module
  co-membership matrix

  **INNOVATION:** The compensatory shift score is the single most novel
  feature in the pipeline. No prior study has simultaneously measured,
  for all RTK pairs, how much one RTK is transcriptionally upregulated in
  response to inhibition of its partner in TNBC tissue. This feature is
  predicted (by SHAP analysis) to be the strongest predictor in the ML
  redundancy model.
  -----------------------------------------------------------------------

  -----------------------------------------------------------------------
  **F1.6 RTK Network Co-Community Function (NCCF-RTK)**

  **MECHANISM:** Constructs the RTK-only PPI subgraph from STRING v12
  (confidence ≥700, experimental + co-expression + text-mining channels).
  Applies Leiden community detection (leidenalg, resolution 0.5--1.5) to
  identify RTK functional communities. Two RTKs in the same community
  share more interaction partners than expected by chance, indicating
  shared pathway membership. Betweenness centrality computed per RTK to
  identify hub RTKs whose loss would most disrupt network flow.

  **INPUTS:** STRING v12 RTK-centred subgraph; NetworkX graph
  construction; leidenalg community detection

  **OUTPUTS:** RTK community partition labels; betweenness/eigenvector
  centrality per RTK; top RTK hub ranking; RTK-RTK co-community matrix
  (binary, 58×58); network figure data for Figure 4

  **INNOVATION:** Community-level analysis reveals RTK pairs that are
  functionally co-embedded in the same signalling module --- a
  systems-level definition of redundancy that is independent of sequence
  or structure. RTK pairs that are both high-KDSI AND co-community are
  the highest-confidence redundancy candidates.
  -----------------------------------------------------------------------

  -----------------------------------------------------------------------
  **F1.7 RTK-RTK Ensemble Redundancy Score (RRS₁)**

  **MECHANISM:** Integrates all five preceding feature sets (F1.1--F1.6)
  into a single ensemble ML redundancy score using Random Forest (500
  trees) + Gradient Boosting (300 estimators) + Elastic Net, weighted
  average. Training labels: \~40 positive pairs (literature-confirmed
  RTK-RTK compensation events) + \~100 negative pairs (RTKs from
  non-overlapping pathways). SHAP TreeExplainer produces per-pair feature
  importance. Output score 0.0--1.0 where 1.0 = maximal RTK-RTK
  redundancy evidence.

  **INPUTS:** Feature matrix from F1.1--F1.6; curated training labels
  from literature; scikit-learn ensemble models; SHAP explainability

  **OUTPUTS:** RRS₁ score for all 1,653 RTK pairs; SHAP waterfall plots
  for top 20 pairs; ranked list top-20 RTK-RTK redundancy pairs; model
  validation metrics (AUC, R², 5-fold CV); Figure 2 heatmap and top-pairs
  summary

  **INNOVATION:** First quantitative, multi-evidence RTK-RTK redundancy
  score for all human RTK pairs. Prior work scores individual pairs ad
  hoc. The SHAP explainability layer means every score has a mechanistic
  justification --- not a black-box prediction.
  -----------------------------------------------------------------------

**RTK-RTK: Key Biological Examples**

  -----------------------------------------------------------------------------------
  **RTK Pair** **Family**     **Known Redundancy       **Expected   **Clinical
                              Mechanism**              RRS₁**       Implication**
  ------------ -------------- ------------------------ ------------ -----------------
  EGFR --      ErbB           ERBB3 lacks kinase       0.88--0.94   Cetuximab
  ERBB3                       activity but forms                    (anti-EGFR)
                              heterodimers with EGFR;               resistance driven
                              ERBB3-PI3K signalling                 by ERBB3
                              fully bypasses EGFR                   upregulation;
                              inhibition via                        co-target with
                              heregulin-driven                      pertuzumab
                              dimerisation                          (anti-ERBB2/3)

  EGFR -- MET  Cross-family   MET amplification drives 0.72--0.82   MET amplification
                              EGFR-independent                      is the dominant
                              PI3K-AKT activation;                  mechanism of EGFR
                              shared substrates GAB1                TKI resistance in
                              and SHC1                              TNBC; EGFR+MET
                                                                    co-inhibition
                                                                    tested in
                                                                    clinical trials

  FGFR1 --     FGFR           Share 70% kinase domain  0.79--0.87   FGFR1
  FGFR2                       identity; both activate               amplification in
                              FRS2-RAS-ERK; ligand                  TNBC; pan-FGFR
                              promiscuity (FGF1/2 bind              inhibitors
                              both)                                 (erdafitinib)
                                                                    target both but
                                                                    resistance
                                                                    emerges via FGFR2
                                                                    upregulation

  MET -- RON   MET family     82% kinase domain        0.83--0.91   RON upregulation
  (MST1R)                     identity; both bind                   documented after
                              HGF-family ligands;                   MET inhibitor
                              activate identical                    treatment;
                              GAB1-PI3K cascade                     highest KDSI RTK
                                                                    pair in the
                                                                    dataset

  AXL -- MERTK TAM family     Share 61% kinase domain  0.74--0.83   AXL inhibitors
                              identity; both bind                   (bemcentinib)
                              GAS6/PROS1; activate                  show MERTK-driven
                              identical SRC-PI3K                    escape; dual
                              signalling; both promote              AXL+MERTK
                              immune evasion in TNBC                inhibition under
                                                                    investigation

  KDR -- FLT1  VEGFR family   Both bind VEGF; share    0.61--0.71   Anti-angiogenic
                              44% kinase domain                     resistance in
                              identity; FLT1 acts as                TNBC partly
                              VEGF decoy but can                    driven by FLT1
                              activate SRC when                     upregulation
                              overexpressed                         after KDR
                                                                    (VEGFR2) blockade
  -----------------------------------------------------------------------------------

+-----------+----------------------------------------------------------+
| **2**     | **NRTK -- NRTK Redundancy**                              |
|           |                                                          |
| HEAD      | *Intra-class compensation between non-receptor tyrosine  |
|           | kinases*                                                 |
+-----------+----------------------------------------------------------+

**Biological Rationale**

Non-receptor tyrosine kinases lack transmembrane domains and reside in
the cytoplasm, nucleus, or inner leaflet of the plasma membrane. They
are activated by RTK signalling, cell adhesion receptors (integrins),
immune receptors (TCR, BCR), and GPCRs. NRTK families (SRC, JAK, ABL,
FAK, ZAP70/SYK, TEC, CSK, FES) compensate for each other through
fundamentally different mechanisms than RTKs:

-   SH2 domain binding promiscuity: NRTKs use SH2 domains to bind
    phosphotyrosine motifs on activated RTKs and scaffolds. Two NRTKs
    with similar SH2 binding specificity (similar pY+1 to pY+4 motif
    preference) compete for the same docking sites and phosphorylate the
    same downstream substrates

-   SH3 domain scaffold co-occupancy: SH3 domains bind proline-rich
    motifs (PxxP). NRTKs with overlapping SH3 ligand selectivity
    co-occupy the same scaffold complexes (e.g. FAK-paxillin-talin focal
    adhesion complex) and can substitute for each other at these nodes

-   Mutual exclusivity in phosphoproteome activation: In some cell
    contexts, SRC and YES1 are mutually exclusive in their activation
    --- when SRC is inhibited, YES1 phosphorylates the identical
    substrate set with no net change in downstream signal. This
    represents the purest form of NRTK redundancy

-   Transactivation through shared adaptors: ABL1 and ABL2 both bind the
    same CRK/CRKL adaptors; FAK and PYK2 both phosphorylate paxillin
    Y118; SRC and FYN both phosphorylate cortactin Y421 --- making these
    pairs functionally interchangeable at the substrate level

**NRTK-NRTK Redundancy Functions**

Seven computational functions quantify NRTK-NRTK redundancy, producing
the NRTK-NRTK Redundancy Score (RRS₂) for all 496 NRTK pairs.

  -----------------------------------------------------------------------
  **F2.1 SH2 Domain Binding Specificity Matrix (SH2-BSM)**

  **MECHANISM:** The most important NRTK-specific function. For each of
  the 32 NRTKs, builds a position weight matrix (PWM) representing the
  pY+1 to pY+4 binding preference of its SH2 domain, derived from
  experimental SPOT peptide array data (where available) or from
  structural contacts in co-crystal structures. KL divergence D₂(i∥∥j)
  between PWMᵢ and PWMⱼ quantifies SH2 specificity similarity. Low KL
  divergence = near-identical substrate recognition = functional
  SH2-level redundancy.

  **INPUTS:** SPOT array data from literature; SH2 domain structural
  contacts (PDB); UniProt SH2 domain boundaries for all 32 NRTKs; MUSCLE
  alignment of SH2 domains

  **OUTPUTS:** 32×32 SH2 KL-divergence matrix; SH2-based redundancy
  clusters (hierarchical clustering of KL matrix); list of NRTK pairs
  with KL divergence \<0.5 (high SH2 specificity overlap); novel
  convergent SH2 pair candidates

  **INNOVATION:** This 32×32 SH2 binding specificity matrix for all human
  NRTKs simultaneously does not exist in the published literature. Prior
  work analyses individual NRTK SH2 domains in isolation. The full matrix
  reveals convergent SH2 pairs --- NRTKs that are phylogenetically
  distant but have evolved identical substrate recognition --- invisible
  to sequence-based analysis.
  -----------------------------------------------------------------------

  -----------------------------------------------------------------------
  **F2.2 SH3 Domain Proline-Rich Motif Overlap Function (SH3-PMO)**

  **MECHANISM:** Extracts SH3 domain sequences for all 32 NRTKs
  (UniProt/InterPro boundaries). Builds PWMs for PxxP motif recognition
  from structural data and known SH3-ligand co-crystal contacts. Pairwise
  SH3 specificity similarity computed as Pearson correlation of PWM
  columns. Two NRTKs with high SH3 similarity co-occupy the same scaffold
  complexes (e.g. FAK--paxillin, SRC--EGFR tail), meaning either can
  occupy the same protein complex and phosphorylate co-localised
  substrates.

  **INPUTS:** SH3 domain sequences (UniProt); co-crystal structure
  SH3-ligand contacts (PDB); known PxxP motif databases

  **OUTPUTS:** 32×32 SH3 similarity matrix; NRTK pairs with overlapping
  SH3 ligand selectivity; scaffold co-occupancy predictions; combined
  SH2+SH3 overlap score per NRTK pair

  **INNOVATION:** SH3 redundancy has been studied pairwise (SRC vs YES1,
  FAK vs PYK2) but never systematically across all 32 human NRTKs. The
  combined SH2+SH3 overlap score is a new composite descriptor of NRTK
  substrate-level redundancy.
  -----------------------------------------------------------------------

  -----------------------------------------------------------------------
  **F2.3 NRTK Shared Phosphosubstrate Overlap Function (SPOF-NRTK)**

  **MECHANISM:** Identical logic to F1.4 but applied to NRTK pairs.
  Computes Jaccard index of phosphosubstrate sets for all 496 NRTK pairs
  using PhosphoSitePlus. Additionally computes the pY site-level overlap:
  for NRTK pair (i,j), counts the number of identical pY sites (same
  protein, same residue) phosphorylated by both. Identical pY site
  phosphorylation is the strongest possible evidence of functional
  redundancy --- it means both NRTKs produce literally the same
  phosphoproteomic output at that site.

  **INPUTS:** PhosphoSitePlus NRTK kinase-substrate-pY data; 32 NRTK
  substrate sets

  **OUTPUTS:** 496-pair Jaccard matrix; identical pY site count per pair;
  NRTK convergence substrates (pY sites phosphorylated by 3+ NRTKs);
  bipartite NRTK-substrate network

  **INNOVATION:** The identical pY site overlap metric is more stringent
  than Jaccard substrate overlap --- it captures true biochemical
  redundancy at single residue resolution, not just protein-level
  substrate sharing. NRTK pairs with high identical pY overlap are the
  most tractable co-targets.
  -----------------------------------------------------------------------

  -----------------------------------------------------------------------
  **F2.4 NRTK Kinase Domain Similarity Function (KDSF-NRTK)**

  **MECHANISM:** HMMER extraction + MUSCLE alignment of kinase domains
  for all 32 NRTKs. Pairwise identity matrix computed. Unlike RTKs where
  extracellular domains drive substrate specificity, for NRTKs the
  SH2/SH3 domains are more determinant --- so kinase domain identity is a
  secondary (but still informative) redundancy predictor. Particularly
  relevant for SRC family NRTKs where kinase domain identity is \>70%
  across members.

  **INPUTS:** 32 NRTK canonical sequences; HMMER PF00069; MUSCLE
  alignment

  **OUTPUTS:** 32×32 NRTK kinase domain identity matrix; SRC family
  sub-clustering; activation loop comparison; gatekeeper residue table
  for NRTKs

  **INNOVATION:** Provides the structural complement to the SH2/SH3
  functional analysis. SRC family NRTKs (SRC, YES1, FYN, LYN, LCK, HCK,
  BLK, FGR) share \>70% kinase domain identity AND near-identical SH2
  specificity --- making them the highest-confidence NRTK redundancy
  group in the dataset.
  -----------------------------------------------------------------------

  -----------------------------------------------------------------------
  **F2.5 NRTK Co-Expression & Compensatory Shift Function (CERF-NRTK)**

  **MECHANISM:** Identical computational logic to F1.5 but applied to 496
  NRTK pairs. Spearman co-expression across TCGA-BRCA TNBC samples.
  Compensatory shift score measuring upregulation of NRTK j when NRTK i
  is pharmacologically inhibited (using GEO dasatinib/saracatinib
  treatment datasets where SRC-family NRTKs are inhibited). Mutual
  exclusivity analysis: tests whether high expression of NRTK i is
  associated with low expression of NRTK j across patients (mutual
  exclusivity is a signature of functional substitution).

  **INPUTS:** TCGA-BRCA TNBC VST counts; GEO SRC-family inhibitor
  datasets; RPPA NRTK activity data

  **OUTPUTS:** 496-pair Spearman matrix; compensatory shift scores;
  mutual exclusivity p-values per pair; NRTK WGCNA module membership

  **INNOVATION:** The mutual exclusivity analysis for NRTK pairs is novel
  in TNBC. Mutual exclusivity of expression implies that the two NRTKs
  perform the same function in different tumours --- a patient-level
  signature of functional redundancy that is directly translatable to
  patient stratification.
  -----------------------------------------------------------------------

  -----------------------------------------------------------------------
  **F2.6 NRTK Network Vulnerability Function (NVF-NRTK)**

  **MECHANISM:** Constructs the NRTK-only PPI subgraph from STRING v12.
  Applies sequential in silico node removal to compute network robustness
  after each NRTK pair is simultaneously removed. The network
  vulnerability score for pair (i,j) = decrease in the size of the
  largest connected component when both i and j are removed
  simultaneously. High vulnerability score = removing this NRTK pair
  maximally disrupts the NRTK signalling network --- identifying the
  highest-priority co-inhibition target pairs.

  **INPUTS:** STRING v12 NRTK subgraph; NetworkX robustness analysis;
  leidenalg community detection

  **OUTPUTS:** Network vulnerability matrix for all 496 NRTK pairs; top
  10 most network-disruptive NRTK pairs; Leiden community labels for NRTK
  modules; NRTK betweenness centrality ranking

  **INNOVATION:** Network vulnerability-based NRTK pair prioritisation
  has not been applied in TNBC. It provides a systems-level rationale for
  combination therapy that is independent of molecular mechanism ---
  identifying NRTK pairs whose co-inhibition would produce the greatest
  downstream signal disruption.
  -----------------------------------------------------------------------

  -----------------------------------------------------------------------
  **F2.7 NRTK-NRTK Ensemble Redundancy Score (RRS₂)**

  **MECHANISM:** Integrates F2.1--F2.6 feature sets into ensemble ML
  score. Same architecture as RRS₁ (Random Forest + Gradient Boosting +
  Elastic Net) but with NRTK-specific training labels (\~25 positive
  pairs: FAK↔PYK2, SRC↔YES1, SRC↔FYN, ABL1↔ABL2, JAK1↔JAK2, etc.) and
  NRTK-specific feature weights. The SH2-BSM (F2.1) KL divergence is
  expected to be the dominant SHAP feature for NRTKs, in contrast to
  compensatory shift which dominates for RTKs.

  **INPUTS:** Feature matrix from F2.1--F2.6; NRTK-specific training
  labels; ensemble ML models; SHAP explainability

  **OUTPUTS:** RRS₂ score for all 496 NRTK pairs; SHAP waterfall for top
  20 NRTK pairs; ranked top-20 NRTK-NRTK redundancy list; Figure 3
  heatmap data

  **INNOVATION:** First quantitative NRTK-NRTK redundancy score
  incorporating SH2 binding specificity as a primary feature. The
  dominance of SH2 similarity (rather than kinase domain identity) in
  SHAP feature importance would constitute a paradigm-shifting finding
  for the NRTK redundancy field.
  -----------------------------------------------------------------------

**NRTK-NRTK: Key Biological Examples**

  -----------------------------------------------------------------------------------
  **NRTK       **Family**   **Known Redundancy        **Expected   **Clinical
  Pair**                    Mechanism**               RRS₂**       Implication**
  ------------ ------------ ------------------------- ------------ ------------------
  SRC -- YES1  SRC family   78% kinase domain         0.91--0.96   Saracatinib (SRC
                            identity; near-identical               inhibitor)
                            SH2 specificity; both                  resistance driven
                            phosphorylate FAK Y925,                by YES1; SRC+YES1
                            paxillin Y118, cortactin               co-inhibition
                            Y421; YES1 upregulated                 required; YES1
                            after SRC inhibition in                amplification in
                            TNBC                                   10--15% TNBC

  FAK -- PYK2  FAK family   66% kinase domain         0.85--0.92   FAK inhibitors
                            identity; both                         (defactinib)
                            phosphorylate paxillin                 activate
                            Y118 and GRB2 Y160; share              PYK2-mediated
                            focal adhesion complex                 survival
                            co-localisation; PYK2                  signalling; dual
                            compensates FAK in                     FAK+PYK2 inhibitor
                            suspension conditions                  (VS-4718) in
                                                                   clinical
                                                                   development

  ABL1 -- ABL2 ABL family   89% kinase domain         0.87--0.93   BCR-ABL imatinib
                            identity; both bind                    resistance partly
                            CRK/CRKL adaptors; ABL2                driven by ABL2;
                            compensates ABL1                       ponatinib targets
                            (BCR-ABL) in imatinib                  both; relevant in
                            resistance                             TNBC with ABL
                                                                   overexpression

  JAK1 -- JAK2 JAK family   38% kinase domain         0.71--0.80   Ruxolitinib
                            identity but nearly                    (JAK1/2 inhibitor)
                            identical substrate                    active in TNBC
                            (STAT3/STAT5                           with IL-6
                            phosphorylation); both                 signalling;
                            activated by shared                    single-JAK
                            cytokine receptors                     inhibition
                            (IL-6R, gp130)                         ineffective due to
                                                                   JAK1↔JAK2
                                                                   compensation

  SRC -- FYN   SRC family   76% kinase domain         0.82--0.89   FYN upregulation
                            identity; both                         in
                            phosphorylate EGFR Y1068               brain-metastatic
                            and FAK Y397; FYN                      TNBC after SRC
                            specifically active in                 inhibition;
                            neural invasion in TNBC                brain-tropic TNBC
                                                                   may require
                                                                   SRC+FYN
                                                                   co-inhibition

  ZAP70 -- SYK ZAP70/SYK    40% kinase domain         0.68--0.77   SYK expressed in
               family       identity; both activated               \~30% TNBC;
                            by immunoreceptor ITAM                 ZAP70-SYK
                            motifs; both                           redundancy limits
                            phosphorylate LAT Y191;                effectiveness of
                            SYK active in basal-like               single-agent SYK
                            TNBC                                   inhibitors
                                                                   (entospletinib) in
                                                                   TNBC
  -----------------------------------------------------------------------------------

+-----------+----------------------------------------------------------+
| **3**     | **RTK -- NRTK Combined Redundancy**                      |
|           |                                                          |
| HEAD      | *Cross-class compensation --- the most clinically        |
|           | dangerous redundancy dimension*                          |
+-----------+----------------------------------------------------------+

**Biological Rationale**

RTK--NRTK cross-class redundancy is mechanistically distinct from
intra-class redundancy and represents the most clinically important
dimension of kinase compensation in TNBC. It arises through four
mechanisms that create independent routes to the same downstream
effectors:

-   RTK→NRTK downstream relay: RTKs activate NRTKs as obligate
    downstream effectors. When the RTK is inhibited, if the NRTK can be
    activated by an alternative upstream signal (another RTK, integrin,
    or GPCR), it sustains the downstream cascade independently. Classic
    example: EGFR inhibition → SRC activated by integrin signalling →
    FAK-PI3K-AKT sustained

-   NRTK transactivation of RTK: NRTKs can phosphorylate RTKs on their
    intracellular activation loop or C-terminal tail, activating the RTK
    independent of its extracellular ligand. SRC phosphorylates EGFR
    Y845 (activation loop), sustaining EGFR kinase activity even when
    the extracellular domain is blocked by cetuximab. This creates a
    circular RTK↔NRTK co-dependency

-   Convergence on shared phosphosubstrates from independent upstream
    routes: The most powerful form of cross-class redundancy. EGFR and
    SRC independently phosphorylate GAB1 Y627 --- blocking one leaves
    the other to maintain GAB1 phosphorylation and PI3K-AKT signalling.
    The substrate phosphorylation state is unchanged despite
    single-kinase inhibition

-   Signalling hub co-occupancy: At multi-protein signalling complexes
    (focal adhesions, endosomes, lipid rafts), RTKs and NRTKs physically
    co-localise and cross-activate each other. The EGFR-SRC-FAK ternary
    complex at focal adhesions is a TNBC-specific signalling hub where
    all three kinases mutually sustain each other\'s activity

**RTK-NRTK Combined Redundancy Functions**

Eight computational functions quantify RTK-NRTK cross-class redundancy,
producing the Combined Redundancy Score (CRS₃) for all 1,856 RTK-NRTK
pairs. This is the most complex and most novel analysis in the pipeline.

  -----------------------------------------------------------------------
  **F3.1 Cross-Class Shared Phosphosubstrate Function (CSPF)**

  **MECHANISM:** The anchor function of Head 3. For each of the 1,856
  RTK-NRTK pairs (i,j) where i ∈ RTK and j ∈ NRTK, computes the Jaccard
  index of their phosphosubstrate sets (PhosphoSitePlus). Additionally
  identifies cross-class convergence substrates --- pY sites
  phosphorylated by both an RTK and an NRTK. These sites represent the
  ultimate nodes of signalling convergence: a cell can maintain their
  phosphorylation regardless of which class of kinase inhibitor is
  applied.

  **INPUTS:** PhosphoSitePlus all kinase-substrate-pY relationships; 58
  RTK substrate sets + 32 NRTK substrate sets

  **OUTPUTS:** 1,856-pair cross-class Jaccard matrix; cross-class
  convergence substrate list (pY sites phosphorylated by ≥1 RTK AND ≥1
  NRTK); RTK-NRTK bipartite substrate network (primary Phase 4 network);
  convergence substrate clinical biomarker candidates

  **INNOVATION:** Cross-class convergence substrates are a new biological
  concept introduced by this pipeline. They represent phosphoproteomic
  nodes that are completely refractory to single-class kinase inhibition.
  Identifying and clinically validating these as resistance biomarkers
  (measurable in liquid biopsy) is a standalone Nature-publishable
  discovery.
  -----------------------------------------------------------------------

  -----------------------------------------------------------------------
  **F3.2 RTK→NRTK Signalling Relay Function (SRF)**

  **MECHANISM:** Constructs directed signalling relay graphs from
  PhosphoSitePlus and literature-curated pathway data. For each RTK,
  identifies which NRTKs are obligate or optional downstream effectors
  (RTK → adaptor → NRTK). Computes relay redundancy: how many independent
  paths exist from RTK i to NRTK j activation, and what is the shortest
  alternative path if the direct relay is blocked. High relay redundancy
  = blocking the RTK does not block the NRTK.

  **INPUTS:** PhosphoSitePlus directed kinase-substrate data;
  KEGG/Reactome RTK pathway topologies; STRING directed interaction data;
  NetworkX shortest-path analysis

  **OUTPUTS:** RTK→NRTK directed relay graph (58×32 adjacency matrix with
  path redundancy weights); top relay pairs with multiple independent
  activation paths; relay bypass score per RTK-NRTK pair

  **INNOVATION:** Quantifying relay path redundancy --- how many
  independent routes connect an RTK to an NRTK --- has not been done
  systematically. High relay redundancy directly predicts that RTK
  inhibition will fail to block the NRTK, explaining why EGFR inhibitors
  don\'t suppress SRC activity in TNBC.
  -----------------------------------------------------------------------

  -----------------------------------------------------------------------
  **F3.3 NRTK→RTK Transactivation Function (TAF)**

  **MECHANISM:** Identifies RTK-NRTK pairs where the NRTK can
  phosphorylate and activate the RTK intracellularly, independent of
  extracellular ligand. Data sources: PhosphoSitePlus (NRTK → RTK pY
  phosphorylation events), literature-curated transactivation
  relationships (SRC→EGFR Y845, SRC→ERBB2 Y877, SRC→MET Y1234/1235).
  Computes a transactivation score based on the number of documented and
  predicted RTK activation loop pY sites that the NRTK can phosphorylate.

  **INPUTS:** PhosphoSitePlus NRTK→RTK phosphorylation events; crystal
  structures showing RTK activation loop accessibility; RTK activation
  loop pY site annotations

  **OUTPUTS:** NRTK→RTK transactivation matrix (32×58); top
  transactivation pairs (SRC→EGFR, SRC→ERBB2, SRC→MET); transactivation
  score per pair; directed edges for relay graph (F3.2)

  **INNOVATION:** Systematic mapping of all NRTK-to-RTK transactivation
  events has not been published. This function reveals which NRTK
  inhibitions would paradoxically reduce RTK activity --- and conversely,
  which NRTKs sustain RTK activity even when extracellular RTK-targeted
  therapies are applied.
  -----------------------------------------------------------------------

  -----------------------------------------------------------------------
  **F3.4 Cross-Class Co-Expression & Compensation Function (CCCF)**

  **MECHANISM:** For all 1,856 RTK-NRTK pairs, computes: (1) Spearman
  co-expression across TCGA-BRCA TNBC samples; (2) cross-class
  compensatory shift score using GEO datasets where RTK inhibitors are
  applied (measuring NRTK upregulation) and NRTK inhibitors are applied
  (measuring RTK upregulation); (3) Lehmann subtype specificity: which
  TNBC subtype (BL1/BL2/M/MSL/IM/LAR) shows the highest RTK-NRTK
  co-expression for each pair, producing a subtype-specific cross-class
  co-expression atlas.

  **INPUTS:** TCGA-BRCA TNBC VST counts; GEO RTK+NRTK inhibitor treatment
  datasets; Lehmann subtype labels from Phase 3

  **OUTPUTS:** 1,856-pair cross-class Spearman matrix; bidirectional
  compensatory shift scores (RTK→NRTK AND NRTK→RTK); subtype-specific
  RTK-NRTK co-expression atlas (6 subtypes × 1,856 pairs); top pairs per
  subtype

  **INNOVATION:** The bidirectional compensatory shift (measuring both
  RTK upregulation after NRTK inhibition AND NRTK upregulation after RTK
  inhibition) is a new metric capturing the circular co-dependency
  between RTK and NRTK classes. The subtype-specific co-expression atlas
  is the first of its kind and directly enables subtype-stratified
  combination therapy design.
  -----------------------------------------------------------------------

  -----------------------------------------------------------------------
  **F3.5 Cross-Class Network Hub & Bridge Function (NHBF)**

  **MECHANISM:** In the full integrated RTK+NRTK+other-protein PPI
  network, identifies \'bridge nodes\' --- kinases (either RTK or NRTK)
  that sit on the shortest paths between the RTK and NRTK communities.
  These bridge kinases are the translators between RTK-driven and
  NRTK-driven signalling. Computes: betweenness centrality for all 90
  kinases in the integrated network; identifies the top bridge RTK-NRTK
  pairs by combined centrality contribution; detects RTK-NRTK pairs that
  are co-members of the same Leiden community in the integrated network.

  **INPUTS:** Full integrated RTK+NRTK PPI network (STRING v12 + PSP);
  NetworkX betweenness centrality; leidenalg on integrated network

  **OUTPUTS:** Bridge kinase ranking (top 20 hub kinases in integrated
  network); RTK-NRTK co-community membership matrix; integrated network
  Leiden partition; Figure 4 network layout data showing RTK-NRTK
  community structure

  **INNOVATION:** Bridge kinases --- those sitting at the interface
  between RTK and NRTK communities --- are the highest-priority
  therapeutic targets in the dataset. They represent nodes whose
  inhibition disrupts both RTK-driven AND NRTK-driven signalling
  simultaneously, overcoming cross-class compensation in a single drug.
  -----------------------------------------------------------------------

  -----------------------------------------------------------------------
  **F3.6 Signalling Hub Co-Occupancy Function (SHCOF)**

  **MECHANISM:** Identifies RTK-NRTK pairs that physically co-localise at
  the same intracellular signalling complex. Data sources: protein
  complex databases (CORUM), co-immunoprecipitation data (BioPlex 3.0),
  proximity ligation assay data from literature. RTK-NRTK pairs
  co-occupying the same complex have direct physical access to each
  other\'s substrates, making cross-class substrate compensation
  structurally inevitable. Computes a co-occupancy score based on the
  number of shared complex memberships.

  **INPUTS:** CORUM protein complex database; BioPlex 3.0 proximity
  interaction data; literature-curated RTK-NRTK co-immunoprecipitation
  events

  **OUTPUTS:** RTK-NRTK co-occupancy matrix (binary, complex-level); top
  co-occupancy pairs; list of shared complexes per pair; co-occupancy
  score as ML feature for F3.8

  **INNOVATION:** Physical co-localisation at the same signalling complex
  is the mechanistic prerequisite for substrate-level cross-class
  compensation. The EGFR-SRC-FAK focal adhesion complex, the
  ERBB2-SRC-PI3K endosomal complex, and the MET-SRC-GAB1 scaffolding
  complex are predicted to emerge as the dominant RTK-NRTK co-occupancy
  hubs in TNBC.
  -----------------------------------------------------------------------

  -----------------------------------------------------------------------
  **F3.7 AI-Grounded Cross-Class Interpretation Function (ACIF)**

  **MECHANISM:** The AI-native function unique to Head 3. For each
  top-ranked RTK-NRTK pair from F3.1--F3.6, the LLM agent (Llama 3.1
  70B + RAG retrieval from Qdrant literature corpus) generates a
  mechanistic interpretation: (1) retrieves top-5 most relevant
  literature passages from the 30,000-vector PubMed corpus; (2)
  synthesises a mechanistic hypothesis explaining the RTK-NRTK
  redundancy; (3) identifies any published evidence contradicting the
  statistical prediction (contradiction = potential novel discovery); (4)
  proposes a testable experimental validation strategy.

  **INPUTS:** Qdrant vector database (30,000+ PubMedBERT embeddings);
  LlamaIndex RAG retriever; LangGraph agent with PSP lookup and KEGG
  pathway tools; statistical results from F3.1--F3.6

  **OUTPUTS:** Per-pair mechanistic interpretation JSON; contradiction
  detection log (predicted redundant pairs where literature says
  \'unrelated\'); novel hypothesis proposals; first-draft manuscript
  Methods/Discussion paragraphs for top RTK-NRTK pairs

  **INNOVATION:** No cancer genomics study has used an AI-native RAG+LLM
  loop to interpret its own statistical results in real time. The
  contradiction detection capability --- finding cases where network
  analysis predicts redundancy but literature says the pair is unrelated
  --- is the most direct route to discovering genuinely new biology.
  -----------------------------------------------------------------------

  -----------------------------------------------------------------------
  **F3.8 RTK-NRTK Combined Redundancy Score (CRS₃)**

  **MECHANISM:** Master integration function. Combines all features from
  F3.1--F3.7 plus the individual RRS₁ and RRS₂ scores into a final
  ensemble ML model specifically trained on cross-class redundancy.
  Training labels: \~50 positive pairs (literature-confirmed RTK-NRTK
  compensatory events: EGFR→SRC, SRC→EGFR transactivation, EGFR-SRC
  shared substrate GAB1, MET-SRC co-activation, etc.) + \~150 negative
  pairs. Produces the master Combined Redundancy Score (CRS₃) for all
  1,856 RTK-NRTK pairs, with full SHAP decomposition showing which
  evidence dimension (substrate overlap? relay path? co-expression? AI
  literature score?) drives each prediction.

  **INPUTS:** All features from F3.1--F3.7; RRS₁ and RRS₂ as additional
  features; cross-class training labels; ensemble ML (Random Forest +
  Gradient Boosting + Elastic Net); SHAP TreeExplainer

  **OUTPUTS:** CRS₃ score for all 1,856 RTK-NRTK pairs; SHAP global
  feature importance (revealing which dimension most predicts cross-class
  redundancy); SHAP waterfall for top 20 RTK-NRTK pairs; master ranked
  list of all 4,005 kinase pairs (RTK-RTK + NRTK-NRTK + RTK-NRTK
  combined); Figure 4 complete redundancy heatmap

  **INNOVATION:** The CRS₃ is the first quantitative, multi-evidence,
  cross-class kinase redundancy score ever computed. By incorporating the
  AI literature grounding score (from F3.7) as a feature alongside
  molecular and network features, it bridges the gap between
  computational prediction and biological validation in a single model.
  The master ranked list of all 4,005 kinase pairs is the definitive
  output of the entire pipeline.
  -----------------------------------------------------------------------

**RTK-NRTK: Key Biological Examples**

  ------------------------------------------------------------------------------------
  **RTK-NRTK   **Redundancy          **Direction**     **Expected   **TNBC Evidence**
  Pair**       Mechanism**                             CRS₃**       
  ------------ --------------------- ----------------- ------------ ------------------
  EGFR -- SRC  Bidirectional: EGFR   RTK↔NRTK          0.92--0.97   EGFR inhibition
               activates SRC via     (circular)                     increases SRC
               GRB2-SOS; SRC                                        activity in TNBC
               transactivates EGFR                                  cell lines; SRC
               Y845; both                                           inhibition
               phosphorylate GAB1                                   increases EGFR
               Y627 and SHC1 Y317                                   Y1068
                                                                    phosphorylation;
                                                                    EGFR+SRC
                                                                    co-inhibition
                                                                    required for
                                                                    sustained response

  EGFR -- FAK  EGFR activates FAK    RTK→NRTK relay    0.78--0.86   FAK Y397
               via SRC relay; FAK                                   phosphorylation
               sustains PI3K-AKT                                    maintained after
               independently; both                                  EGFR inhibition in
               phosphorylate                                        TNBC; gefitinib +
               paxillin                                             defactinib
                                                                    combination active
                                                                    in EGFR-low TNBC

  MET -- SRC   MET activates SRC via RTK↔NRTK          0.83--0.90   MET inhibitor
               GAB1; SRC can         (circular)                     (crizotinib)
               activate MET                                         resistance driven
               Y1234/Y1235                                          by SRC activation;
               (activation loop                                     MET+SRC
               transactivation);                                    co-inhibition
               share GAB1/GRB2                                      synergistic in
               substrates                                           MET-amplified TNBC

  FGFR1 -- SRC FGFR1 activates SRC   RTK→NRTK relay    0.74--0.83   FGFR inhibitor
               via FRS2-GRB2; SRC                                   (erdafitinib)
               independently                                        resistance in TNBC
               sustains ERK via RAS;                                associated with
               converge on SHC1 Y317                                SRC
                                                                    hyperactivation;
                                                                    FGFR1+SRC
                                                                    co-inhibition
                                                                    proposed

  AXL -- SRC   AXL activates SRC via RTK→NRTK relay    0.76--0.84   AXL inhibition in
               GAS6; SRC sustains                                   TNBC triggers
               AKT via PI3K                                         SRC-STAT3
               independently; both                                  compensatory
               phosphorylate STAT3                                  activation;
               Y705                                                 AXL+SRC combined
                                                                    inhibition shows
                                                                    synergy in
                                                                    AXL-high TNBC

  ERBB2 -- SRC SRC phosphorylates    NRTK→RTK          0.81--0.89   SRC inhibition
               ERBB2 Y877            transactivation                reduces ERBB2
               (activation loop);                                   phosphorylation in
               ERBB2 activates SRC                                  TNBC with ERBB2
               via GRB7; circular                                   low expression;
               co-dependency at                                     SRC+ERBB2
               ERBB2-enriched                                       co-targeting
               vesicles                                             active even in
                                                                    HER2-low TNBC

  MET -- FAK   MET activates FAK via RTK↔NRTK scaffold 0.79--0.87   MET-FAK
               SRC relay and direct                                 co-activation
               phosphorylation; FAK                                 drives invasion in
               sustains MET                                         mesenchymal TNBC
               expression via NF-κB;                                subtype; focal
               converge on paxillin                                 adhesion complex
               complex                                              is the
                                                                    co-occupancy hub
  ------------------------------------------------------------------------------------

**Integrated Scoring: How All Three Heads Combine**

The three redundancy scores (RRS₁, RRS₂, CRS₃) are produced
independently but feed into a master ranked list covering all 4,005
kinase pair combinations. The scores are not averaged --- they are
maintained as separate dimensions to preserve the mechanistic
specificity of each head. A pair can score high on one head and low on
another, and the SHAP decomposition for each score reveals exactly why.

  -------------------------------------------------------------------------------------------
  **Score**     **Pairs       **Features Used**        **Expected Top      **Clinical
                Scored**                               Pair**              Output**
  ------------- ------------- ------------------------ ------------------- ------------------
  RRS₁          1,653 RTK     KDSI, gatekeeper,        MET -- RON (82%     Top 5 RTK-RTK
  (RTK-RTK)     pairs         structure RMSD, pocket   KDSI + compensatory pairs proposed as
                              Tanimoto, substrate      upregulation +      dual-RTK
                              Jaccard, co-expression,  shared GAB1         inhibition targets
                              compensatory shift,      substrates)         in TNBC
                              network community                            

  RRS₂          496 NRTK      SH2 KL-divergence, SH3   SRC -- YES1         Top 5 NRTK-NRTK
  (NRTK-NRTK)   pairs         similarity, pY site      (near-identical SH2 pairs proposed as
                              overlap, kinase domain   specificity +       dual-NRTK
                              identity, co-expression, confirmed           inhibition
                              mutual exclusivity,      compensatory        targets; YES1
                              network vulnerability    upregulation +      inhibitor
                                                       mutual exclusivity  development
                                                       across TNBC         prioritised
                                                       patients)           

  CRS₃          1,856         All F3.1--F3.7           EGFR -- SRC         Top 5 RTK-NRTK
  (RTK-NRTK)    cross-class   features + RRS₁/RRS₂ as  (bidirectional      pairs proposed as
                pairs         meta-features + AI       transactivation +   cross-class
                              literature score         highest shared      co-inhibition
                                                       substrate count +   targets;
                                                       strongest AI        convergence
                                                       literature          substrates
                                                       support + highest   proposed as
                                                       compensatory shift  phosphoproteomic
                                                       in both directions) resistance
                                                                           biomarkers

  Master Ranked All 4,005     Union of all three score To be determined by Single definitive
  List          pairs         dimensions + head label  data --- the novel  ranked list
                                                       discovery will be   submitted as
                                                       the top-ranked pair supplementary
                                                       that has NOT been   data; top 10 pairs
                                                       previously          become primary
                                                       described as        manuscript
                                                       redundant           findings
  -------------------------------------------------------------------------------------------

  -------------- ------------------------------------------------------------
  **THE MASTER   The three-head framework produces three independently
  INNOVATION**   validated redundancy scores that together cover the entire
                 kinase pair universe. The cross-head analysis ---
                 identifying kinase pairs that score highly on ALL THREE
                 heads simultaneously (i.e. an RTK and NRTK that also share
                 substrates with other RTKs and NRTKs) --- identifies
                 multi-connected hub pairs whose inhibition would
                 simultaneously disrupt RTK-RTK, NRTK-NRTK, and RTK-NRTK
                 signalling. These represent the ultimate TNBC co-targets: a
                 single combination therapy addressing all three redundancy
                 dimensions at once.

  -------------- ------------------------------------------------------------

**Complete Function Reference**

All 22 redundancy functions implemented across the three analytical
heads:

  -----------------------------------------------------------------------------------------------
  **Function   **Name**           **Head**    **Phase**   **Key Output**     **Novel
  ID**                                                                       Contribution**
  ------------ ------------------ ----------- ----------- ------------------ --------------------
  F1.1         Kinase Domain      RTK-RTK     2           58×58 KDSI matrix  First complete RTK
               Sequence Identity                                             KDSI matrix

  F1.2         Gatekeeper Residue RTK-RTK     2           Gatekeeper matrix; All-RTK gatekeeper
               Identity                                   inhibitor          comparison
                                                          cross-reactivity   unpublished

  F1.3         Structural Binding RTK-RTK     2           RMSD + pocket      All-RTK simultaneous
               Pocket Similarity                          Tanimoto matrices  structural
                                                                             comparison

  F1.4         Shared             RTK-RTK     1+4         Jaccard matrix;    RTK convergence
               Phosphosubstrate                           convergence        substrate concept
               Overlap (RTK)                              substrates         

  F1.5         Co-Expression &    RTK-RTK     3           Spearman + shift   Bidirectional
               Compensatory Shift                         matrix             compensatory shift
               (RTK)                                                         metric

  F1.6         Network            RTK-RTK     4           Community labels;  RTK community-level
               Co-Community (RTK)                         betweenness        redundancy
                                                          ranking            definition

  F1.7         RTK-RTK Ensemble   RTK-RTK     6           RRS₁ for 1,653     First quantitative
               Score (RRS₁)                               pairs + SHAP       RTK-RTK redundancy
                                                                             score

  F2.1         SH2 Binding        NRTK-NRTK   2           32×32 SH2          Entire NRTK SH2
               Specificity Matrix                         KL-divergence      specificity matrix
                                                          matrix             unpublished

  F2.2         SH3 Proline-Rich   NRTK-NRTK   2           SH3 similarity     Systematic SH3
               Motif Overlap                              matrix             comparison across
                                                                             all NRTKs

  F2.3         Shared             NRTK-NRTK   1+4         Jaccard +          Identical pY site
               Phosphosubstrate                           identical pY site  metric
               Overlap (NRTK)                             matrix             (single-residue
                                                                             resolution)

  F2.4         NRTK Kinase Domain NRTK-NRTK   2           32×32 NRTK KDSI    SRC family
               Similarity                                 matrix             cross-redundancy
                                                                             quantification

  F2.5         NRTK Co-Expression NRTK-NRTK   3           Spearman + mutual  NRTK mutual
               & Mutual                                   exclusivity matrix exclusivity as
               Exclusivity                                                   redundancy signature

  F2.6         NRTK Network       NRTK-NRTK   4           Vulnerability      Network-robustness
               Vulnerability                              matrix; top        based NRTK pair
                                                          co-inhibition      prioritisation
                                                          pairs              

  F2.7         NRTK-NRTK Ensemble NRTK-NRTK   6           RRS₂ for 496       First quantitative
               Score (RRS₂)                               pairs + SHAP       NRTK-NRTK redundancy
                                                                             score

  F3.1         Cross-Class Shared RTK-NRTK    1+4         Cross-class        Cross-class
               Phosphosubstrate                           convergence        convergence
                                                          substrates         substrate concept
                                                                             (new)

  F3.2         RTK→NRTK           RTK-NRTK    4           Relay path         Systematic relay
               Signalling Relay                           redundancy scores  path quantification
                                                                             unpublished

  F3.3         NRTK→RTK           RTK-NRTK    1+4         Transactivation    All-NRTK to all-RTK
               Transactivation                            matrix             transactivation map

  F3.4         Cross-Class        RTK-NRTK    3           Bidirectional      Subtype-specific
               Co-Expression &                            shift + subtype    RTK-NRTK
               Compensation                               atlas              co-expression atlas

  F3.5         Cross-Class        RTK-NRTK    4           Bridge kinase      Bridge kinase
               Network Hub &                              ranking            concept for
               Bridge                                                        cross-class
                                                                             targeting

  F3.6         Signalling Hub     RTK-NRTK    1+4         Co-occupancy       Systematic RTK-NRTK
               Co-Occupancy                               matrix; shared     complex co-occupancy
                                                          complex list       map

  F3.7         AI-Grounded        RTK-NRTK    5           Mechanistic        First AI-RAG
               Cross-Class                                paragraphs;        interpretation of
               Interpretation                             contradiction log  its own results

  F3.8         RTK-NRTK Combined  RTK-NRTK    6           CRS₃ for 1,856     First cross-class
               Score (CRS₃)                               pairs; master      kinase redundancy
                                                          ranked list        score ever
  -----------------------------------------------------------------------------------------------

*RTK--NRTK Redundancy Functions Reference · February 2026*
