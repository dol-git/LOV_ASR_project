## Overview

This repository documents an undergraduate research project on the evolutionary design principles of LOV-domain photoreceptor proteins.

The project combines ancestral sequence reconstruction (ASR), phylogenetic analysis, indel-aware sequence refinement, AlphaFold-based structural evaluation, and literature-based construct design to select ancestral LOV-domain candidates for future biochemical and spectroscopic experiments.

The central aim is not simply to reconstruct ancestral sequences, but to understand how a conserved LOV light-sensing core diversified into different signaling systems.

The main biological questions are:

1. How did LOV-domain photochemistry evolve?
2. How did dark recovery kinetics diversify among LOV proteins?
3. How did different output modules become linked to a common LOV sensor core?
4. Can evolutionary principles of LOV proteins inform synthetic biology and optogenetic protein design?

This repository currently reflects the computational design and construct selection stage. Experimental validation will be added in future updates.

---

## Background

LOV domains are blue-light-sensing protein modules that bind flavin cofactors such as FMN. Upon blue-light irradiation, a conserved cysteine residue forms a covalent adduct with the flavin cofactor, triggering structural changes that are transmitted to output domains.

Although the LOV core is conserved, different LOV-containing proteins use this core in distinct ways.

Representative examples include:

| Protein family | Representative protein | Architecture | Signaling strategy |
|---|---|---|---|
| Plant phototropin-like LOV2 | AsLOV2 | LOV2 + Jα / kinase-associated region | Monomeric conformational change |
| LOV-STAS | YtvA | LOV + STAS | Dimeric signaling system |
| LOV-HTH | EL222 | LOV + HTH DNA-binding domain | Light-induced dimerization and DNA binding |

This project focuses on comparing these major LOV architectures from an evolutionary perspective.

---

## Research Concept

The guiding idea of this project is:

> A conserved LOV sensory core can be evolutionarily coupled to different output mechanisms.

Therefore, the project asks whether ancestral LOV-domain sequences can reveal common design principles behind:

- photochemistry
- dark recovery kinetics
- structural switching
- oligomeric-state changes
- output-domain coupling

The selected ancestral constructs are intended for future comparison with modern LOV proteins by spectroscopy and structural analysis.

---

## Repository Status

Current milestone:

```text
Computational construct design and selection completed.

This repository currently contains:

selected ancestral construct candidates
tree files used for phylogenetic comparison
sequence identity summaries
construct design outputs
literature-based construct boundary rationale
documentation of the selection strategy

Experimental results are not yet included.

Future updates may include:

expression and purification results
UV-visible absorption spectra
dark recovery kinetics
comparison with modern LOV proteins
oligomeric-state analysis
final interpretation of ancestral functional properties
Project Workflow

The workflow was organized as follows:

Collection of LOV-containing protein sequences
Domain architecture classification
Dataset filtering and separation into plant-like and bacterial-like groups
Multiple sequence alignment using several MAFFT strategies
Phylogenetic tree inference using IQ-TREE
Tree topology comparison across alignment methods
Ancestral sequence reconstruction using ConsistASR
Indel-aware refinement of reconstructed ancestral sequences
Candidate node extraction
Sequence convergence and pairwise identity analysis
AlphaFold-based structural evaluation
Literature-based construct boundary mapping
Final ancestral construct selection

The overall strategy was to avoid relying on a single alignment, tree, or reconstructed sequence. Instead, candidate robustness was evaluated across multiple computational conditions.

Dataset Construction

Representative LOV proteins were used as starting points:

AsLOV2: plant phototropin-like LOV2 representative
YtvA: bacterial LOV-STAS representative
EL222: bacterial LOV-HTH representative

The project initially considered full-length LOV-containing proteins rather than only isolated LOV domains. This was important because the research question concerns not only the LOV core itself, but also how different output modules became associated with the LOV sensor.

Sequences were classified into major domain architecture groups, including:

phototropin-like LOV / kinase-associated proteins
YtvA-like LOV-STAS proteins
EL222-like LOV-HTH proteins
bacterial LOV-HK-like proteins

After filtering, separate plant-like and bacterial-like datasets were prepared for phylogenetic analysis and ASR.

Domain Architecture Classification

Domain architecture classification was used to avoid mixing biologically different LOV systems too casually.

The major categories were:

Category	Description
Phototropin-like LOV	Plant-like LOV proteins associated with kinase regions
LOV-STAS	Bacterial YtvA-like proteins containing LOV and STAS domains
LOV-HTH	EL222-like proteins containing LOV and helix-turn-helix DNA-binding domains
LOV-HK	LOV proteins associated with histidine kinase-like regions

This classification was important because different output domains imply different evolutionary and functional constraints.

Alignment Strategy

Multiple sequence alignment was performed using several MAFFT strategies:

MAFFT AUTO
MAFFT L-INS-i
MAFFT E-INS-i
MAFFT G-INS-i

The use of multiple alignment strategies was necessary because LOV-containing proteins often include variable linker regions and different output domains. These regions can strongly affect alignment quality and downstream ASR results.

Candidate reliability was assessed by checking whether similar ancestral nodes, sequences, and structures were recovered across different alignment conditions.

In this project, alignment robustness was treated as an important part of candidate confidence.

Phylogenetic Analysis

Phylogenetic trees were inferred using IQ-TREE with model selection.

Separate plant-like and bacterial-like trees were constructed. Tree topology comparisons were then used to evaluate whether major clades were stable across alignment strategies.

The trees/ directory contains representative tree files and topology comparison outputs.

Tree comparison was especially important because ASR depends directly on tree topology. If a candidate ancestral node only appears under one unstable alignment condition, it is less reliable as an experimental target.

Ancestral Sequence Reconstruction

Ancestral sequence reconstruction was performed using ConsistASR.

After ASR, indel-aware refinement was applied to account for insertion and deletion events. This step produced both gap-containing and gap-stripped ancestral sequences.

Indel-aware reconstruction was important because construct design requires realistic protein sequences, not only residue substitutions on a fixed alignment length.

The final candidate selection did not rely on only one reconstructed sequence. Instead, candidates were evaluated across multiple alignment and reconstruction conditions.

Candidate Selection Strategy

Final candidates were selected using the following criteria:

Sequence convergence across alignment strategies
Structural consistency in AlphaFold predictions
Conservation of the LOV core and photoactive cysteine region
Agreement with experimentally validated modern LOV construct boundaries
Interpretability as representative ancestral nodes
Suitability for experimental expression and purification

Consensus sequences were tested, but they were not used as the final strategy.

Instead, representative best-node sequences were selected. This was because best-node sequences preserve clearer phylogenetic meaning, while consensus sequences can become artificial constructs that do not correspond to a specific reconstructed ancestral state.

Sequence Identity and Convergence Analysis

Candidate sequences were compared within each family using pairwise identity calculations.

The goal was not to force all sequences to become identical, but to identify candidates that were central or representative among multiple reconstruction conditions.

Pairwise identity summaries are stored in:

constructs/*_pairwise.tsv
constructs/*_pairwise.summary.tsv

These analyses helped support the choice of final representative nodes.

Structural Evaluation

AlphaFold-based structural prediction was used to evaluate whether candidate sequences produced plausible LOV-domain folds.

The main structural checks were:

preservation of the LOV core fold
structural convergence across alignment-derived candidates
absence of obvious structural collapse
reasonable behavior of linker or terminal regions
consistency with known modern LOV-domain structures

For some groups, especially LOV-HK-like candidates, predicted structures were harder to interpret. In contrast, the HTH, STAS, and plant LOV2-like candidates showed clearer structural behavior and were therefore prioritized.

Why Consensus Sequences Were Not Used

Consensus sequences were initially generated for major groups, but they were not selected as final experimental candidates.

The reasons were:

Consensus sequences do not correspond to a specific ancestral node.
Ambiguous regions can hide real evolutionary differences.
Best-node sequences preserve phylogenetic interpretability.
Structural predictions supported the use of representative nodes.
The final goal is evolutionary interpretation, not only sequence optimization.

Therefore, the final constructs were selected from specific ancestral nodes rather than from group-level consensus sequences.

Final Selected Constructs

Three ancestral LOV-domain constructs were selected for experimental consideration.

Construct	Source group	Reference boundary	Length
hth_Node83_EL222_1_144	EL222-like / LOV-HTH	EL222 residues 1–144	122 aa
plant_Node62_AsLOV2_404_560	AsLOV2-like / plant LOV2	AsLOV2 residues 404–560	157 aa
stas_Node227_YtvA_20_147	YtvA-like / LOV-STAS	YtvA residues 20–147	128 aa

The final construct sequences are stored under:

constructs/construct_design/outputs/

A combined final FASTA file may also be stored as:

constructs/final_LOV_construct_candidates.fasta
Final Construct Sequences
EL222-like / LOV-HTH candidate
>hth_Node83_EL222_1_144 mapped_from_EL222_ref_1_144 length=122
MAAGDADESLSAMIANSPIAAVITNPRLPDNPIVACNEAFCDLTGYSRDEIIGRNCRFLSGPDTEPELTEQIREAIRERRPVLVEILNYKKDGTPFRNAVMVAPIFDEDGELEYFLGSQVEV
AsLOV2-like / plant LOV2 candidate
>plant_Node62_AsLOV2_404_560 mapped_from_AsLOV2_ref_404_560 length=157
LATTLERIEKNFVITDPRLPDNPIIFASDSFLELTEYSREEILGRNCRFLQGPETDRATVRKIRDAIDNQRDVTVQLINYTKSGKKFWNLFHLQPMRDQKGEVQYFIGVQLDGSEHVEEDTAKEGAKLVKETADNVDEAVRELPDANLKPEDLWANH
YtvA-like / LOV-STAS candidate
>stas_Node227_YtvA_20_147 mapped_from_YtvA_ref_20_147 length=128
LDHTRVAVVITDPSQPDNPIIYANKGFTELTGYSEDEVIGRNCRFLQGEDTDPETVDKIREAIREKEPVTVELLNYKKDGTPFWNELHIDPIYIEDKYYFVGTQKDITEQKEAEQKLEDKHKEIERLS
Literature-Based Construct Boundaries

Construct boundaries were chosen by mapping ancestral candidates onto experimentally used modern LOV constructs.

This was done to avoid arbitrary trimming and to make the selected constructs experimentally realistic.

EL222-like construct

Reference protein:

EL222 from Erythrobacter litoralis HTCC2594

Reported experimental constructs:

Full-length EL222: residues 1–222
EL222 LOV fragment: residues 1–144

Expression system reported in literature:

Vector: pET-28a
Host: Escherichia coli Rosetta2(DE3)

Reference:

Photoinduced dimerization of a photosensory DNA-binding protein EL222 and its LOV domain

The ancestral HTH candidate was mapped to the EL222 LOV fragment boundary. Because the ancestral sequence did not contain the same N-terminal extension as modern EL222, the final mapped construct became 122 aa rather than 144 aa. The conserved LOV core and photoactive cysteine region were retained.

AsLOV2-like construct

Reference protein:

AsLOV2 from Avena sativa phototropin

Reference boundary:

AsLOV2 residues 404–560

Expression system reported in literature:

Vector: pET-19b
Host: Escherichia coli BL21(DE3)

Reference:

Reduction midpoint potential of a paradigm light–oxygen–voltage receptor and its modulation by methionine residues

The selected plant-like ancestral construct includes the LOV2 core and surrounding region corresponding to the experimentally used AsLOV2 region.

YtvA-like construct

Reference protein:

YtvA from Bacillus subtilis

Reported experimental construct:

isolated YtvA LOV domain, approximately residues 1–127 in the reference study

Selected reference mapping:

YtvA residues 20–147

Expression system reported in literature:

Vector: pET-41a
Host: Escherichia coli CmpX13

Reference:

Time-Resolved X-Ray Solution Scattering Reveals the Structural Photoactivation of a Light-Oxygen-Voltage Photoreceptor

The selected STAS-like ancestral construct was trimmed around the YtvA LOV-domain region and preserves the conserved LOV core.

Why LOV-HK Was Not Selected

LOV-HK-like proteins were examined during candidate selection, but they were not selected for the first experimental construct set.

The main reasons were:

Domain boundaries were harder to define confidently.
Linker and kinase-associated regions showed greater variability.
Structural predictions were less straightforward to interpret.
Some candidates showed lower confidence in non-LOV regions.
Literature-supported construct boundaries were less clear for the selected ancestral candidates.

This does not mean that LOV-HK proteins are unimportant. Rather, they likely require a separate and more careful analysis before experimental construct design.

For the first experimental set, HTH, STAS, and plant LOV2-like candidates were prioritized because they showed clearer sequence, structural, and literature-based support.

Functional Diversity of Selected LOV Systems

Although all selected constructs share a conserved LOV core, they represent distinct signaling strategies.

Group	Modern representative	Signaling strategy
AsLOV2-like	AsLOV2	monomeric conformational change, including Jα helix rearrangement
YtvA-like	YtvA	dimeric LOV-STAS signaling
EL222-like	EL222	light-induced dimerization and DNA binding

This diversity is central to the project. The selected constructs are intended to help test how different signaling mechanisms evolved from a common LOV-domain framework.

Repository Structure
.
├── README.md
├── constructs
│   ├── bacterial
│   ├── plant
│   ├── consensus
│   ├── construct_design
│   ├── metadata
│   ├── *_pairwise.tsv
│   └── *_pairwise.summary.tsv
├── data
├── docs
├── figures
├── results
├── scripts
│   └── lov_architecture_classifier.py
└── trees
    ├── *_v2.treefile
    ├── bacterial_only
    ├── plant_only
    └── comparison
Directory Notes
constructs/

Contains selected candidate sequences, intermediate candidate FASTA files, consensus tests, pairwise identity summaries, and final construct design outputs.

Important files include:

constructs/construct_design/outputs/
constructs/*_pairwise.summary.tsv
constructs/*_pairwise.tsv
trees/

Contains representative tree files and topology comparison outputs.

Important contents include:

trees/*_v2.treefile
trees/comparison/
trees/bacterial_only/
trees/plant_only/
scripts/

Contains selected scripts that were useful during the project.

This repository is not currently intended to be a fully automated software package. Some scripts were used as temporary analysis utilities and may require cleanup before reuse.

figures/

Reserved for summary figures such as:

phylogenetic tree comparison
selected candidate structures
construct boundary diagrams
sequence alignment snapshots
docs/

Reserved for extended notes, methodology descriptions, and literature summaries.

Reproducibility Notes

This repository is currently a result-oriented research record.

It is not yet a fully reproducible pipeline from raw sequence collection to final construct selection.

The most important reproducible elements are:

selected construct sequences
tree files
pairwise identity summaries
construct boundary mappings
literature-based expression system information
documentation of selection criteria

Future versions may include:

cleaned command-line scripts
example datasets
standardized environment files
reproducible workflow scripts
automatic figure generation
Current Limitations

This project is still in progress.

Current limitations include:

Experimental validation has not yet been completed.
Some analysis scripts were written as temporary utilities.
Oligomeric-state analysis is planned but not yet finalized.
LOV-HK-like candidates require additional analysis.
Spectroscopic characterization remains future work.

These limitations will be addressed in later stages of the project.

Planned Experimental Direction

The selected constructs are intended for future expression and purification.

Planned experimental analyses include:

protein expression screening
purification of ancestral LOV constructs
UV-visible absorption spectroscopy
blue-light response measurement
dark recovery kinetics
comparison with modern LOV proteins
possible oligomeric-state analysis

The long-term goal is to connect sequence evolution, structural behavior, photochemistry, and signaling mechanism.

Version History
v1: Construct design and selection stage

Current version.

Includes:

candidate selection strategy
final ancestral construct sequences
phylogenetic tree files
pairwise sequence identity summaries
construct boundary rationale
literature-supported expression system information
Future v2: Experimental validation stage

Planned update after experiments.

Expected additions:

expression results
purification results
spectroscopy data
dark recovery kinetics
comparison with modern LOV proteins
updated figures and interpretation
Summary

This project selected three representative ancestral LOV-domain constructs corresponding to major LOV signaling architectures:

EL222-like LOV-HTH
AsLOV2-like plant LOV2
YtvA-like LOV-STAS

The candidates were selected through a multi-step process combining phylogenetic robustness, ASR, indel-aware refinement, structural prediction, sequence comparison, and literature-based construct boundary mapping.

The next major step is experimental validation.