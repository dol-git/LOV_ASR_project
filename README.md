# LOV ASR Project

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

LOV domains are blue-light-sensing protein modules that bind flavin cofactors such as FMN.

Upon blue-light irradiation, a conserved cysteine residue forms a covalent adduct with the flavin cofactor, triggering structural changes that are transmitted to output domains.

Although the LOV core is conserved, different LOV-containing proteins use this core in distinct ways. In this project, three representative architecture classes were used as reference systems for construct design.

## Representative architectures

| Protein family | Representative protein | Architecture | Signaling strategy |
|---|---|---|---|
| Plant phototropin-like LOV2 | AsLOV2 | LOV2 + Jalpha / kinase-associated region | Monomeric conformational change |
| LOV-STAS | YtvA | LOV + STAS | LOV-STAS linker-mediated signaling |
| LOV-HTH | EL222 | LOV + HTH DNA-binding domain | Light-regulated DNA binding |

---

## Research Concept

A conserved LOV sensory core can be evolutionarily coupled to different output mechanisms.

This project uses ancestral sequence reconstruction to compare how LOV-domain proteins may have diversified across different signaling architectures.

---

## Repository Status

Current milestone: computational construct design and selection completed.

## This repository currently contains

- selected ancestral construct candidates
- tree files used for phylogenetic comparison
- sequence identity summaries
- construct design outputs
- literature-based construct boundary rationale
- documentation of the selection strategy

## Not yet included

- experimental validation
- spectroscopy
- dark recovery kinetics

---

## Project Workflow

1. Collection of LOV-containing protein sequences
2. Domain architecture classification
3. Dataset filtering
4. Multiple sequence alignment using MAFFT
5. Phylogenetic tree inference using IQ-TREE
6. Tree topology comparison
7. Ancestral sequence reconstruction using ConsistASR
8. Indel-aware sequence refinement
9. Candidate extraction
10. Sequence convergence analysis
11. AlphaFold-based structural evaluation
12. Literature-based construct boundary mapping
13. Final expression construct selection

---

## Repository Structure

| Directory | Contents |
|---|---|
| `constructs/` | Candidate ancestral sequences, consensus sequences, pairwise identity tables, and construct-design working files |
| `constructs/construct_design/` | Working files used for reference-based construct extraction |
| `results/construct_design/` | Clean construct-design records prepared for GitHub documentation |
| `results/construct_design/final_order_fastas/` | Final FASTA records for selected WT and ancestral constructs |
| `trees/` | IQ-TREE tree files and topology comparison outputs |
| `scripts/` | Utility scripts for classification and analysis |
| `figures/` | Figures for summaries and future reports |
| `docs/` | Additional documentation |
| `data/` | Data placeholder directory |

---

## Dataset Construction

Representative LOV proteins used as reference systems:

- AsLOV2
- YtvA
- EL222

Sequences were classified into architecture groups including:

- phototropin-like LOV
- LOV-STAS
- LOV-HTH

---

## Domain Architecture Classification

| Category | Description |
|---|---|
| Phototropin-like LOV | Plant-like LOV proteins related to phototropin LOV2 systems |
| LOV-STAS | YtvA-like proteins containing LOV and STAS modules |
| LOV-HTH | EL222-like proteins containing LOV and HTH DNA-binding modules |

---

## Alignment Strategy

Multiple alignment strategies were used to evaluate robustness:

- MAFFT AUTO
- MAFFT L-INS-i
- MAFFT E-INS-i
- MAFFT G-INS-i

The goal was to identify ancestral candidates that were not artifacts of a single alignment condition.

---

## Phylogenetic Analysis

Phylogenetic trees were inferred using IQ-TREE.

Tree topology comparisons were used to evaluate whether candidate ancestral nodes were stable across alignment strategies.

Tree files and comparison outputs are stored in:

- `trees/`
- `trees/comparison/`

---

## Ancestral Sequence Reconstruction

Ancestral sequence reconstruction was performed using ConsistASR.

Indel-aware refinement was used to avoid artificially overextended ancestral sequences and to generate more realistic candidate proteins for structure prediction and experimental design.

---

## Candidate Selection Strategy

Candidate ancestral sequences were evaluated using the following criteria:

- sequence convergence across alignment conditions
- structural consistency
- preservation of the conserved LOV cysteine
- compatibility with known reference proteins
- literature-supported construct boundaries
- phylogenetic interpretability
- suitability for future expression and spectroscopy

---

## Construct Design Strategy

The construct design step was not limited to extracting the conserved LOV core.

For each architecture, the construct boundary was chosen to include the relevant C-terminal helix or linker region that may participate in signal transmission.

The final design used both:

1. Literature-supported reference boundaries
2. The actual retained ancestral residues after reference-coordinate extraction

This was important because ancestral sequences can contain gaps relative to modern reference proteins. Therefore, the length of an ancestral construct does not always match the length of the corresponding wild-type reference construct.

Detailed construct-design records are stored in:

- `results/construct_design/README.md`
- `results/construct_design/final_order_fastas/`

The working construct-design files are stored in:

- `constructs/construct_design/`

---

## Final Selected Construct Set

| Construct name | Reference boundary | Length | Status |
|---|---:|---:|---|
| `WT_AsLOV2_404_560` | AsLOV2 404-560 | 157 aa | Ordered |
| `AncPlant_Node62_AsLOV2_404_560eq` | AsLOV2 404-560 equivalent | 157 aa | Ordered |
| `WT_YtvA_20_147` | YtvA 20-147 | 128 aa | Prepared |
| `AncSTAS_Node227_YtvA_20_147eq` | YtvA 20-147 equivalent | 128 aa | Prepared |
| `WT_EL222_1_163` | EL222 1-163 | 163 aa | Prepared |
| `AncHTH_Node83_EL222_1_170eq` | EL222 1-170 equivalent | 148 aa | Prepared |

---

## Construct Boundary Rationale

### AsLOV2 / plant-like LOV2

The wild-type AsLOV2 construct was defined as residues 404-560.

This region contains the LOV2 core and the C-terminal Jalpha-containing region. The ancestral plant-like construct was extracted from Node62 using the same AsLOV2 404-560 reference coordinates.

Both the wild-type and ancestral AsLOV2-like constructs have been submitted for ordering.

### YtvA / LOV-STAS

The YtvA construct was defined as residues 20-147.

This region contains the LOV core and the downstream LOV-STAS linker / Jalpha-like segment while excluding the STAS domain body.

The ancestral Node227 sequence preserved the corresponding linker region with minimal gaps, so the YtvA 20-147-equivalent construct was retained.

### EL222 / LOV-HTH

The EL222 reference boundary was defined as:

- LOV domain: residues 1-144
- connector Jalpha-helix: residues 145-163
- HTH domain: residues 164-222

Therefore, the wild-type EL222 construct was defined as residues 1-163.

For the ancestral HTH candidate, strict extraction using the EL222 1-163 boundary retained a relatively short C-terminal connector segment due to alignment gaps. Therefore, EL222 1-163, 1-170, and 1-175 equivalent ancestral constructs were compared.

The EL222 1-170-equivalent ancestral construct was selected because it preserved the Jalpha-like charged/helical connector segment more completely while minimizing extension into the HTH domain.

---

## Final FASTA Records

Final construct FASTA files are stored in:

- `results/construct_design/final_order_fastas/all_WT_and_ancestral_order_constructs.fasta`
- `results/construct_design/final_order_fastas/order_YtvA_EL222_WT_Ancestor_constructs_renamed.fasta`

---

## Sequence Identity Analysis

Sequence identity summaries are stored in:

- `constructs/hth_pairwise.tsv`
- `constructs/hth_pairwise.summary.tsv`
- `constructs/plant_pairwise.tsv`
- `constructs/plant_pairwise.summary.tsv`
- `constructs/stas_pairwise.tsv`
- `constructs/stas_pairwise.summary.tsv`

---

## Structural Evaluation

AlphaFold-based structural evaluation was used to assess:

- fold conservation
- structural stability
- domain integrity
- C-terminal helix or linker plausibility

Further structure-prediction inputs for the final construct set will be prepared after final ordering decisions.

---

## Summary

The project has selected three ancestral LOV-domain candidates and their corresponding wild-type reference constructs:

- AsLOV2-like plant LOV2 candidate
- YtvA-like LOV-STAS candidate
- EL222-like LOV-HTH candidate

The AsLOV2 wild-type and ancestral constructs have already been submitted for ordering. The YtvA and EL222 wild-type/ancestral construct sets are prepared as the next ordering candidates.

---

## Next Steps

1. Submit the YtvA and EL222 WT/ancestor constructs for ordering.
2. Read key papers on AsLOV2, YtvA, and EL222 signal-transmission helices/linkers.
3. Prepare AlphaFold or ColabFold inputs for the final construct set.
4. After synthesis, proceed to cloning, expression screening, purification, and spectroscopy.
