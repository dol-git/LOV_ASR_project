# Software Environment

This file records software used in the computational workflow. Versions are
reported only when confirmed from local commands. Entries without a confirmed
version are marked **to be confirmed**.

## Verified Local Environment

The following versions were queried from the local `ConsistASR` Conda
environment:

| Tool | Version | Role |
|---|---|---|
| Python | 3.12.3 | Script execution and workflow support |
| Biopython | 1.86 | FASTA, alignment, and sequence processing |
| MAFFT | 7.525 | Multiple sequence alignment |
| IQ-TREE | 2.0.7 | Phylogenetic inference and ancestral reconstruction |
| RAxML-NG | 1.2.2 | Phylogenetic/indel-aware workflow support |
| RAxML | 8.2.12 | Indel-aware workflow support |

The locally available HMMER command reports:

| Tool | Version | Role |
|---|---|---|
| HMMER (`hmmscan`) | 3.4 | Pfam domain annotation |

## Workflow Components

| Component | Version | Notes |
|---|---|---|
| ConsistASR | to be confirmed | Indel-aware ASR workflow developed in the lab and used to refine reconstructed ancestral sequences |
| AlphaFold | to be confirmed | Structure-prediction screening used during candidate evaluation |
| Pfam database | to be confirmed | Profile-HMM database used with HMMER for domain annotation |

## Reproducibility Notes

- The verified versions above describe the currently available local
  environment and should be checked against archived run logs before a formal
  release.
- Exact commands, input checksums, database release identifiers, and
  environment export files remain to be consolidated.
- AlphaFold predictions are computational screening results, not experimental
  structure validation.

## HMR analysis snapshot (reviewed 2026-10-21)

These pins are transcribed from the existing HMR analysis requirements file,
not a newly installed or freshly tested environment. Existing independent
validation metadata also records NumPy, MDAnalysis and MDTraj versions.
Do not substitute the ASR environment above for this analysis environment.

| Package/application | Recorded version | Use |
|---|---|---|
| NumPy | 2.4.6 | Coordinates, RMSF, contact summaries |
| MDAnalysis | 2.10.0 | PSF/DCD/PDB I/O and distance arrays |
| MDTraj | 1.11.2 | Eight-state DSSP |
| SciPy | 1.17.1 | Pairwise distances for medoid/FPS |
| Matplotlib | 3.11.0 | Historical plots (not regenerated or bundled) |
| ChimeraX | 1.11.1 | Native viewer; historical headless validation |
| Python used for HMR | to be confirmed | ASR Python version is not evidence for HMR |
| GENESIS executable used for production | to be confirmed | Exact run binary/build not certified here |

The original script prepends its local `software/` directory for MDTraj;
that vendor directory is excluded. A future separate workspace must supply the
recorded packages and external inputs. No dependency installation or version
compatibility run was performed for this update.
