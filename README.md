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

Although the LOV core is conserved, different LOV-containing proteins use this core in distinct ways.

### Representative architectures

| Protein family | Representative protein | Architecture | Signaling strategy |
|---|---|---|---|
| Plant phototropin-like LOV2 | AsLOV2 | LOV2 + Jα / kinase-associated region | Monomeric conformational change |
| LOV-STAS | YtvA | LOV + STAS | Dimeric signaling system |
| LOV-HTH | EL222 | LOV + HTH DNA-binding domain | Light-induced dimerization and DNA binding |

---

## Research Concept

> A conserved LOV sensory core can be evolutionarily coupled to different output mechanisms.

---

## Repository Status

**Current milestone:**  
Computational construct design and selection completed.

### This repository currently contains

- selected ancestral construct candidates  
- tree files used for phylogenetic comparison  
- sequence identity summaries  
- construct design outputs  
- literature-based construct boundary rationale  
- documentation of the selection strategy  

### Not yet included

- experimental validation  
- spectroscopy  
- kinetics  

---

## Project Workflow

1. Collection of LOV-containing protein sequences  
2. Domain architecture classification  
3. Dataset filtering  
4. Multiple sequence alignment (MAFFT)  
5. Phylogenetic tree inference (IQ-TREE)  
6. Tree topology comparison  
7. Ancestral sequence reconstruction (ConsistASR)  
8. Indel-aware refinement  
9. Candidate extraction  
10. Sequence convergence analysis  
11. AlphaFold evaluation  
12. Construct boundary mapping  
13. Final selection  

---

## Dataset Construction

Representative LOV proteins used:

- AsLOV2  
- YtvA  
- EL222  

Sequences were classified into:

- phototropin-like LOV  
- LOV-STAS  
- LOV-HTH  
- LOV-HK  

---

## Domain Architecture Classification

| Category | Description |
|---|---|
| Phototropin-like LOV | Plant-like LOV proteins |
| LOV-STAS | YtvA-like proteins |
| LOV-HTH | EL222-like proteins |
| LOV-HK | Histidine kinase-associated LOV proteins |

---

## Alignment Strategy

- MAFFT AUTO  
- MAFFT L-INS-i  
- MAFFT E-INS-i  
- MAFFT G-INS-i  

Multiple alignment strategies were used to ensure robustness.

---

## Phylogenetic Analysis

Trees were inferred using IQ-TREE.

The `trees/` directory contains tree files and comparison outputs.

---

## Ancestral Sequence Reconstruction

Performed using ConsistASR.

Indel-aware refinement was applied to generate realistic sequences.

---

## Candidate Selection Strategy

Criteria:

- sequence convergence  
- structural consistency  
- conserved cysteine  
- literature agreement  
- phylogenetic interpretability  

---

## Sequence Identity Analysis

Stored in:


constructs/_pairwise.tsv
constructs/_pairwise.summary.tsv


---

## Structural Evaluation

AlphaFold was used to evaluate:

- fold conservation  
- structural stability  
- domain integrity  

---

## Final Selected Constructs

| Construct | Length |
|---|---|
| hth_Node83_EL222_1_144 | 122 aa |
| plant_Node62_AsLOV2_404_560 | 157 aa |
| stas_Node227_YtvA_20_147 | 128 aa |

---

## Final Construct Sequences

### EL222-like


MAAGDADESLSAMIANSPIAAVITNPRLPDNPIVACNEAFCDLTGYSRDEIIGRNCRFLSGPDTEPELTEQIREAIRERRPVLVEILNYKKDGTPFRNAVMVAPIFDEDGELEYFLGSQVEV


### AsLOV2-like


LATTLERIEKNFVITDPRLPDNPIIFASDSFLELTEYSREEILGRNCRFLQGPETDRATVRKIRDAIDNQRDVTVQLINYTKSGKKFWNLFHLQPMRDQKGEVQYFIGVQLDGSEHVEEDTAKEGAKLVKETADNVDEAVRELPDANLKPEDLWANH


### YtvA-like


LDHTRVAVVITDPSQPDNPIIYANKGFTELTGYSEDEVIGRNCRFLQGEDTDPETVDKIREAIREKEPVTVELLNYKKDGTPFWNELHIDPIYIEDKYYFVGTQKDITEQKEAEQKLEDKHKEIERLS


---

## Structural Comparison

![Structure](figures/structure_comparison.png)

---

## Summary

Three ancestral LOV constructs were selected:

- EL222-like  
- AsLOV2-like  
- YtvA-like  

---

## Next Step

Experimental validation (expression, spectroscopy, kinetics)
