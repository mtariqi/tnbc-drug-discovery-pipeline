**PRESENTER SCRIPTS**

*Decoding Kinase Redundancy & Crosstalk in TNBC*

Slide-by-Slide Speaker Notes

Md Tariqul Islam · UAMS · April 2026

**How to use this document**

Each section below corresponds to one slide in the presentation. The
script is written in full sentences to be read aloud or used as the
basis for natural delivery. Suggested durations and presenter tips are
included to help you pace and prepare for likely questions.

**Total estimated presentation time: 25--35 minutes (excluding Q&A).**

**SLIDE 1 Title Slide**

*⏱ Suggested time: 1--2 minutes*

Good \[morning/afternoon\], Professor \[Name\], and thank you for the
opportunity to present this doctoral research proposal.

My name is Md Tariqul Islam, and I am working under the supervision of
Dr. Sayem Miah in the Department of Biochemistry and Molecular Biology
at the University of Arkansas for Medical Sciences.

The title of this proposal is: Decoding Kinase Redundancy and Crosstalk
in Triple-Negative Breast Cancer --- A Comprehensive Multi-Omics
Bioinformatics and Data Science Framework for RTK and NRTK Target
Discovery.

In plain terms, this project asks a simple but clinically critical
question: why do cancer drugs that target individual kinase proteins
consistently fail in triple-negative breast cancer --- and what would a
rational, data-driven approach to combination therapy look like?

This is a 24-month, fully computational research programme drawing on 11
gigabytes of publicly available multi-omics data from The Cancer Genome
Atlas, specifically the TCGA-BRCA cohort.

I will walk you through the scientific rationale, our hypotheses, the
computational methods, and the expected scientific contributions of this
work.

  -----------------------------------------------------------------------
  **💡 PRESENTER TIP:** *Speak slowly and make eye contact. Let the title
  settle before you begin --- don\'t rush into content.*

  -----------------------------------------------------------------------

**SLIDE 2 The Clinical Problem**

*⏱ Suggested time: 2--3 minutes*

Let me start with the clinical landscape so we are all grounded in the
problem this research is trying to solve.

Triple-negative breast cancer --- TNBC --- is defined by the absence of
estrogen receptor, progesterone receptor, and HER2 expression. That
definition is not just biological --- it is also a therapeutic death
sentence for many patients, because it forecloses every available
targeted therapy we currently have in the breast cancer toolkit.

TNBC accounts for 15 to 20 percent of all breast cancer diagnoses
globally, but it is responsible for a disproportionate share of cancer
deaths. The median overall survival at Stage IV is below 18 months. This
is a disease that kills quickly.

It disproportionately affects younger women and Black women ---
communities that already face systemic barriers to healthcare. This is
not just a scientific problem; it is a public health equity problem.

EGFR is overexpressed in 50 to 70 percent of TNBC tumours. FGFR1 and
FGFR2 are amplified in up to 22 percent of cases. AXL, MET, and DDR
kinases drive invasion and immune evasion. By all molecular logic, these
should be excellent drug targets. And yet --- approved targeted kinase
inhibitors have zero clinical track record of sustained benefit as
monotherapy in TNBC.

The reason is kinase redundancy. When you block one kinase, the tumour
simply re-routes its signalling through a related kinase or through a
cross-class partner. The downstream oncogenic output is maintained. The
drug appears to work in isolation, but the cancer doesn\'t care.

Systematic quantification of these redundancy relationships is precisely
what this programme will deliver --- and it is the precondition for
rational combination therapy design.

  -----------------------------------------------------------------------
  **💡 PRESENTER TIP:** *Pause after the survival statistic --- it is
  stark and your audience needs a moment to absorb it. The last bullet is
  your transition into the science, so land it with emphasis.*

  -----------------------------------------------------------------------

**SLIDE 3 Research Hypotheses**

*⏱ Suggested time: 2 minutes*

This proposal is structured around four testable hypotheses, each
corresponding to a distinct dimension of kinase redundancy.

Hypothesis 1 addresses RTK-to-RTK redundancy: we propose that subsets of
the 58 receptor tyrosine kinases expressed in TNBC form functionally
redundant transcriptional modules --- groups where inhibiting any single
member is insufficient to suppress the shared downstream signalling
output.

Hypothesis 2 addresses NRTK-to-NRTK redundancy: we propose that the 32
non-receptor tyrosine kinases active in TNBC are organised into
hierarchical regulatory clusters whose intra-cluster functional overlap
directly underpins adaptive resistance to monotherapy.

Hypothesis 3 is the most clinically dangerous scenario --- cross-class
RTK to NRTK crosstalk: we propose that specific RTK-NRTK pairs maintain
direct, independent correlations in TNBC --- not just because they share
upstream activators, but because they actively communicate and
compensate for each other.

And the Translational Hypothesis, which is the ultimate goal of the
entire programme: that by integrating network centrality, machine
learning essentiality scores, survival data, and druggability evidence,
we can identify a small, tractable set of RTK-NRTK co-target pairs that
are simultaneously network-essential, clinically prognostic, and
confirmed CRISPR-essential in TNBC cell lines.

These are the pairs that would be most likely to yield durable clinical
responses in combination kinase inhibitor trials.

  -----------------------------------------------------------------------
  **💡 PRESENTER TIP:** *Read each hypothesis slowly and deliberately.
  These are the scientific stakes --- not just boxes to check. The
  Translational Hypothesis is the one with direct clinical impact, so
  give it extra weight.*

  -----------------------------------------------------------------------

**SLIDE 4 The Three-Head Redundancy Framework**

*⏱ Suggested time: 3 minutes*

The analytical architecture of this programme is built around what I
call the Three-Head Redundancy Framework --- three distinct levels at
which kinases compensate for each other in TNBC, each requiring its own
computational logic.

Head 1 covers RTK-to-RTK redundancy. There are 58 receptor tyrosine
kinases in the human kinome, which gives us 1,653 unique pairwise
combinations to evaluate. RTKs compensate through trans-phosphorylation
--- when you block EGFR, for example, ERBB3 can form alternate
heterodimers that sustain the same downstream cascade. They also share
ligands, share adapter substrates like GAB1 and SHC1, and trigger
compensatory transcriptional upregulation. The combined redundancy score
for all RTK pairs is called RRS₁.

Head 2 covers NRTK-to-NRTK redundancy --- 496 pairwise combinations
among the 32 non-receptor tyrosine kinases. NRTKs are cytoplasmic
enzymes that compensate through fundamentally different mechanisms: SH2
domain binding promiscuity --- where two NRTKs recognise the same
phosphotyrosine motifs --- and SH3 scaffold co-occupancy, where they
physically substitute for each other at the same multi-protein
complexes. The score here is RRS₂.

Head 3 is cross-class RTK-to-NRTK redundancy --- 1,856 pairs --- and
this is the dimension that is most dangerous clinically and most novel
scientifically. RTKs and NRTKs can sustain each other\'s signalling
through completely independent upstream routes. The canonical example is
EGFR and SRC: block EGFR with cetuximab, and SRC --- activated via
integrin signalling --- independently sustains AKT and ERK through
shared phosphosubstrates. These two kinases don\'t even need to talk to
each other directly. The score for cross-class pairs is CRS₃.

The total analytical scope of this framework is just over 4,000 kinase
pairs, all scored, ranked, and SHAP-explained.

  -----------------------------------------------------------------------
  **💡 PRESENTER TIP:** *Point to each column as you discuss it. The pair
  counts --- 1,653, 496, and 1,856 --- are memorable numbers that give
  the audience a sense of the scale of the analysis. Emphasise Head 3 as
  the most novel and most dangerous.*

  -----------------------------------------------------------------------

**SLIDE 5 Data Acquisition & Multi-Layer QC**

*⏱ Suggested time: 2--3 minutes*

Before any biological inference is drawn, this programme invests
substantial effort in data quality --- an aspect that is frequently
underemphasised in published multi-omics studies but is foundational to
the validity of every downstream conclusion.

All data are publicly available through the NCI Genomic Data Commons
under the TCGA-BRCA project. The dataset is approximately 11 gigabytes
and spans six data modalities: RNA-seq gene-level counts using the
STAR-HTSeq workflow aligned to GRCh38; FPKM and FPKM-UQ normalised
expression matrices; somatic mutation calls in MAF format from three
independent callers; gene-level copy number variation calls from
GISTIC2; full clinical and biospecimen metadata; and Illumina HM450K DNA
methylation beta-value matrices.

TNBC samples are curated using intersection criteria: ER-negative,
PR-negative, and HER2-negative by IHC and FISH. We anticipate 155 to 200
confirmed TNBC specimens after curation.

The QC framework has seven layers, each targeting a distinct failure
mode. Layer 1 checks raw alignment quality --- mapping rates, duplicate
fractions, insert size distributions. Layer 2 assesses gene-level count
integrity and library depth. Layer 3 uses RSeQC Transcript Integrity
Numbers to flag RNA degradation. Layer 4 implements DESeq2 VST
normalisation, ComBat-seq batch correction, and PVCA variance
partitioning. Layer 5 applies consensus somatic variant filtering and
CNV-expression concordance checks. Layer 6 uses PCA, UMAP, Mahalanobis
distance, and hierarchical clustering for multi-modal outlier detection.
And Layer 7 performs RNA-DNA concordance checks to detect potential
sample swaps.

Every sample receives a Composite Sample Quality Score --- a single
quantitative metric that enables auditable, principled inclusion
decisions. This score is itself a methodological contribution of the
programme.

  -----------------------------------------------------------------------
  **💡 PRESENTER TIP:** *Professors appreciate rigour. Spend time on the
  QC slide --- it signals that this is serious, reproducible science, not
  just a bioinformatics fishing expedition. The Composite Sample Quality
  Score is a novel contribution in itself.*

  -----------------------------------------------------------------------

**SLIDE 6 RTK--RTK Redundancy Functions (Head 1)**

*⏱ Suggested time: 3 minutes*

Head 1 implements seven computational functions that together produce
the RTK-RTK Redundancy Score --- RRS₁ --- for all 1,653 RTK pairs.

F1.1 is the Kinase Domain Sequence Identity function. Using HMMER to
extract kinase domains from all 58 RTKs and MUSCLE for multiple sequence
alignment, we compute a 58-by-58 pairwise identity matrix. High identity
indicates that both kinases have near-identical catalytic machinery and
will likely phosphorylate the same substrate tyrosine motifs. This is
the first time this matrix has been computed for all 58 human RTKs
simultaneously.

F1.2 extracts the gatekeeper residue --- the amino acid that controls
access to the ATP-binding back pocket --- for all 58 RTKs. Two kinases
sharing a gatekeeper are cross-inhibited by the same clinical compounds,
which means a drug targeting one may inadvertently leave the other
partially active.

F1.3 performs 3D structural superposition of kinase domain binding
pockets using ProDy and fpocket. Pairs with low structural RMSD and high
Tanimoto similarity in pocket descriptors are structurally redundant ---
the same inhibitor scaffold will bind both.

F1.4 computes the Jaccard index of phosphosubstrate sets from
PhosphoSitePlus for all RTK pairs, and importantly identifies
convergence substrates --- pY sites phosphorylated by three or more RTKs
simultaneously. These sites cannot be blocked by inhibiting any single
RTK. We anticipate that GAB1, IRS1, SHC1, and GRB2 will emerge as the
most clinically critical convergence substrates.

F1.5 is our most novel single feature. It combines Spearman
co-expression across TCGA-BRCA TNBC samples with a Compensatory Shift
Score --- a metric that directly measures, for every RTK pair, how much
RTK-j is transcriptionally upregulated when RTK-i is pharmacologically
suppressed, using publicly available drug-resistance GEO datasets. SHAP
analysis predicts this will be the strongest predictor in the ensemble
model.

F1.6 uses Leiden community detection on the STRING protein interaction
network to identify RTK pairs that are functionally co-embedded in the
same signalling module --- a systems-level definition of redundancy
independent of sequence or structure.

F1.7 integrates all five feature sets using a Random Forest, Gradient
Boosting, and Elastic Net ensemble, with SHAP explainability providing a
mechanistic justification for every score. The output is a ranked list
of the top 20 highest-confidence RTK-RTK redundancy pairs.

  -----------------------------------------------------------------------
  **💡 PRESENTER TIP:** *You don\'t need to explain every function in
  equal depth --- walk through F1.1 to F1.4 at a steady pace, then slow
  down for F1.5 (the compensatory shift score) because that is your most
  novel feature, and finish strongly on F1.7 which is the integration
  step.*

  -----------------------------------------------------------------------

**SLIDE 7 NRTK--NRTK Redundancy Functions (Head 2)**

*⏱ Suggested time: 3 minutes*

Head 2 applies seven analogous but mechanistically distinct functions to
the 496 NRTK pairs, producing the RRS₂ score.

The anchor function is F2.1 --- the SH2 Domain Binding Specificity
Matrix. For NRTKs, it is not primarily the kinase domain that determines
substrate specificity --- it is the SH2 domain, which binds
phosphotyrosine motifs on activated receptors and scaffolds. F2.1 builds
a position weight matrix for the SH2 pY+1 to pY+4 binding preference of
each NRTK from experimental SPOT array data and structural co-crystal
contacts, then computes KL divergence between every pair of PWMs. Low KL
divergence means near-identical substrate recognition --- and therefore
functional SH2-level redundancy.

This 32-by-32 SH2 binding specificity matrix for all human NRTKs
simultaneously does not exist in the published literature. More
importantly, it can reveal convergent SH2 pairs --- NRTKs that are
phylogenetically distant but have evolved identical substrate
recognition --- which are completely invisible to sequence-based
analysis.

F2.2 extends this to SH3 domain promiscuity: NRTKs with overlapping SH3
ligand selectivity can physically substitute for each other at the same
multi-protein scaffolds like the FAK-paxillin-talin focal adhesion
complex.

F2.3 applies phosphosubstrate overlap, now at single pY site resolution
--- the most stringent possible evidence of biochemical redundancy. If
SRC and YES1 both phosphorylate FAK at tyrosine 925, they are
functionally interchangeable at that exact molecular address.

F2.5 adds a mutual exclusivity analysis across TNBC patients: if high
expression of NRTK-i is statistically associated with low expression of
NRTK-j across the cohort, this is a patient-level signature of
functional substitution. This is directly translatable to patient
stratification in clinical trials.

F2.6 computes a Network Vulnerability Score by sequential in silico
removal of NRTK pairs from the STRING interaction network. Pairs whose
co-removal maximally disrupts network connectivity are the
highest-priority co-inhibition candidates --- a systems rationale for
combination therapy independent of molecular mechanism.

The key biological examples --- SRC and YES1 with an expected RRS₂ of
0.91 to 0.96, FAK and PYK2 with 0.85 to 0.92, and ABL1 and ABL2 with
0.87 to 0.93 --- are all clinically validated redundancy pairs that
serve as positive training labels for the ensemble model.

  -----------------------------------------------------------------------
  **💡 PRESENTER TIP:** *Make sure to emphasise F2.1 as your most
  innovative contribution in Head 2. The SH2 story --- that
  phylogenetically distant NRTKs can have identical substrate recognition
  --- is counterintuitive and scientifically exciting. Let that land.*

  -----------------------------------------------------------------------

**SLIDE 8 RTK--NRTK Cross-Class Redundancy (Head 3)**

*⏱ Suggested time: 2--3 minutes*

Head 3 is the most clinically dangerous and most scientifically novel
dimension of this framework --- cross-class redundancy between RTKs and
NRTKs.

Unlike intra-class redundancy, RTK-NRTK compensation does not require
sequence similarity or even shared binding pockets. These two classes of
kinases can independently sustain the same oncogenic output through
fundamentally different upstream routes.

The first mechanism is the downstream relay. When an RTK is inhibited,
an NRTK that was originally downstream of that RTK can be re-activated
by an entirely different upstream signal --- a growth factor receptor,
an integrin, or a GPCR --- and sustain the downstream cascade
independently. This is why blocking EGFR does not stop SRC activity: SRC
has multiple activators.

The second mechanism is transactivation --- and this one is particularly
insidious. SRC can phosphorylate EGFR at tyrosine 845, the activation
loop residue, sustaining EGFR kinase activity even when the
extracellular domain is fully occupied by a blocking antibody like
cetuximab. The drug physically occupies the receptor but the receptor
remains active from the inside.

The third mechanism is shared phosphosubstrate convergence --- the most
powerful form. EGFR and SRC independently phosphorylate GAB1 at tyrosine
627. If you block EGFR, SRC continues to phosphorylate GAB1. The
phosphorylation state of GAB1 is unchanged. PI3K-AKT signalling
continues. From the cancer\'s perspective, nothing happened.

The fourth mechanism is signalling hub co-occupancy. At focal adhesions,
endosomes, and lipid rafts, RTKs and NRTKs physically co-localise in
ternary complexes. The EGFR-SRC-FAK complex in TNBC is the textbook
example: all three kinases mutually sustain each other\'s activity.
Blocking one shifts activity to the others.

Head 3 will compute the CRS₃ score for all 1,856 RTK-NRTK pairs using
eight computational functions --- including a cross-class shared
phosphosubstrate function, an NRTK-RTK transactivation network, and a
cross-class network bridge analysis to identify which RTK-NRTK pairs
serve as critical bridges between the two signalling layers.

  -----------------------------------------------------------------------
  **💡 PRESENTER TIP:** *Head 3 is conceptually the hardest but
  clinically the most compelling. Use the EGFR-SRC example specifically
  --- it is concrete, clinically familiar to your audience, and perfectly
  illustrates why this problem is so dangerous.*

  -----------------------------------------------------------------------

**SLIDE 9 Integrated Analytical Pipeline**

*⏱ Suggested time: 2 minutes*

The programme runs for 24 months across three phases, each with a
distinct scientific function.

Phase I, spanning Months 1 through 6, is entirely devoted to data
science and exploratory data analysis. This includes the seven-layer QC
framework I described, DESeq2 variance-stabilising transformation and
ComBat-seq batch correction, and comprehensive univariate and
multivariate profiling of all 90 kinase genes using PCA, UMAP, and
t-SNE. No biological hypotheses are tested in Phase I. The goal is to
arrive at a curated, high-quality data matrix that can support robust
inference.

Phase II, Months 7 through 18, is the core redundancy analysis. This is
where we implement all three heads of the framework: WGCNA for RTK and
NRTK co-expression modules, ARACNE and VIPER for mutual information
network inference, Graphical Lasso partial correlation for the bipartite
RTK-NRTK crosstalk network, and the full RRS₁, RRS₂, and CRS₃ ensemble
scoring with SHAP explainability.

Phase III, Months 19 through 24, is external validation and
translational output. Top-ranked target pairs are cross-validated in
independent GEO cohorts and CCLE transcriptomic data, and CRISPR
essentiality evidence is integrated from DepMap. Cox proportional
hazards survival modelling links redundancy scores to patient outcomes.
The final deliverable is a fully documented, reproducible Snakemake
pipeline and a suite of open-access Jupyter and R notebooks released to
the research community.

  -----------------------------------------------------------------------
  **💡 PRESENTER TIP:** *This is your methodology overview slide ---
  don\'t over-explain each bullet. Your goal is to convey that this is a
  rigorous, sequential, reproducible programme. Emphasise that Phase I is
  dedicated entirely to data science before any biological conclusions
  are drawn --- professors often challenge this point in multi-omics
  proposals.*

  -----------------------------------------------------------------------

**SLIDE 10 Expected Outcomes & Scientific Impact**

*⏱ Suggested time: 2 minutes*

Let me close by summarising what this programme will deliver and why it
matters.

The primary deliverable is a ranked, druggability-filtered list of
RTK-NRTK co-target pairs for combination therapy development --- the
first such list derived from a systematic, multi-evidence computational
framework applied to all human RTK and NRTK pairs simultaneously.

Supporting that are four datasets with no published precedent: the first
58-by-58 RTK redundancy matrix with SHAP-explained scores; the first
32-by-32 NRTK redundancy matrix driven by SH2 binding specificity; the
first 1,856-pair cross-class CRS₃ scored network; and compensatory shift
scores for every kinase pair, measuring real-world transcriptional
compensation in response to pharmacological inhibition.

From a scientific novelty standpoint, there are four claims this
programme is positioned to make. First, it delivers the first
quantitative, multi-evidence redundancy scores for all human RTK and
NRTK pairs simultaneously --- replacing ad hoc pairwise studies with a
systematic atlas. Second, the Compensatory Shift Score is a genuinely
new biomarker concept --- a signal extracted from drug-resistance
transcriptomics that directly measures adaptive re-routing in TNBC
tissue. Third, if SHAP analysis confirms that SH2 domain similarity is a
stronger predictor of NRTK redundancy than kinase domain identity, this
would constitute a paradigm shift in how we understand NRTK substrate
selection. And fourth, the cross-class CRS₃ framework defines
cross-class redundancy as a computable, rankable property --- which it
has never been before.

All code, pipelines, and notebooks will be released open-access. This
programme does not just produce results --- it produces community
infrastructure.

Thank you. I am happy to take any questions.

  -----------------------------------------------------------------------
  **💡 PRESENTER TIP:** *End with confidence. You have covered a lot of
  ground --- this slide is your chance to synthesise it. Speak to the
  novelty claims directly, especially the compensatory shift score and
  the SH2 specificity matrix. These are the claims that will draw
  reviewer scrutiny and also reviewer excitement.*

  -----------------------------------------------------------------------
