# LOV ASR Project

## Overview

This repository contains scripts and workflows for the analysis of LOV-domain photoreceptor proteins, with a focus on phylogenetic classification and ancestral sequence reconstruction (ASR).

---

## Objective

The goal of this project is to investigate the evolutionary design principles of LOV-domain proteins by:

* collecting LOV-family protein sequences
* classifying proteins based on domain architecture
* preparing datasets for phylogenetic analysis
* reconstructing ancestral LOV-domain sequences
* comparing extant and ancestral proteins in downstream analyses

---

## Repository Structure

```text
LOV_ASR_project/
├── scripts/   # Python and shell scripts for analysis
├── data/      # Input sequence files and metadata
├── results/   # Output files from classification and analysis
├── figures/   # Plots and images
├── docs/      # Notes and report-related files
└── README.md
```

---

## Methods

Protein classification is performed using HMMER domain annotation (Pfam) followed by rule-based classification.

Each protein is categorized based on:

* presence of PAS domain (required for LOV-like classification)
* presence of output domains (kinase, STAS, HTH, etc.)
* presence of unwanted redox-related domains

---

## Usage

```bash
python scripts/lov_architecture_classifier.py \
    --domtblout data/example.domtblout \
    --fasta data/example.fasta \
    --outprefix results/lov_arch
```

---

## Input

* `domtblout`: HMMER output file (`--domtblout`)
* `fasta`: protein sequences corresponding to the search

---

## Output

### Summary table

* `*.summary.tsv`
* Contains classification results for all proteins

### FASTA outputs

* `*.keep.fasta`: high-confidence LOV proteins
* `*.review.fasta`: candidates requiring manual inspection
* `*.exclude_ids.txt`: filtered proteins

### Subtype-specific outputs

* `lov_stas`, `lov_hk_bact`, `lov_hth`, etc.

---

## Workflow

```bash
# Step 1: Domain annotation
hmmscan --domtblout results/domtblout.txt Pfam-A.hmm data/sequences.fasta

# Step 2: Classification
python scripts/lov_architecture_classifier.py \
    --domtblout results/domtblout.txt \
    --fasta data/sequences.fasta \
    --outprefix results/lov_arch

# Step 3: Inspect results
less results/lov_arch.summary.tsv
```

---

## Notes

* Classification is currently based on domain presence only
* Domain order and coverage are not yet considered (planned improvement)
* This repository represents an early-stage pipeline for ASR-based protein analysis
