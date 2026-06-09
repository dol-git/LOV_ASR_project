# Project Status

Last documentation review: June 9, 2026

## Current Phase

The repository represents the completed computational candidate-selection and
construct-design phase of the LOV-domain ASR graduation research project.
The broader project connects this computational work to experimental comparison
of ancestral and extant LOV proteins. Laboratory work is ongoing and
preliminary.

## High-Level Timeline

- Early April 2026: computational environment setup, LOV sequence collection,
  initial alignment, phylogenetic analysis, and ASR workflow testing.
- Mid to late April 2026: dataset refinement, comparison of four MAFFT
  strategies, IQ-TREE/ASR analysis, indel-aware refinement, and
  AlphaFold-based candidate screening.
- Late April to early May 2026: cross-alignment convergence analysis,
  literature-guided construct boundary selection, and finalization of the
  current ancestral/extant construct set.
- Current experimental phase: preliminary expression, purification, and
  spectroscopy work is ongoing for the AsLOV2 ancestral/extant pair, with
  repeated measurements and optimization in progress.

## Completed Computational Work

- LOV-containing sequence collection and architecture-based filtering
- separation of plant phototropin-like and bacterial LOV datasets
- MAFFT AUTO, L-INS-i, E-INS-i, and G-INS-i alignment comparisons
- IQ-TREE phylogenetic inference
- pairwise topology comparison across alignment strategies
- ancestral sequence reconstruction and indel-aware refinement
- extraction of candidate plant, LOV-STAS, LOV-HTH, and LOV-HK ancestors
- sequence convergence analysis using pairwise identity
- AlphaFold-based computational screening
- reference-guided construct-boundary mapping
- preparation of six final extant/ancestral construct designs

## Selected Designs

The current computationally selected set contains:

- AsLOV2 extant and plant-like ancestral Node62 constructs
- YtvA extant and LOV-STAS ancestral Node227 constructs
- EL222 extant and LOV-HTH ancestral Node83 constructs

The AsLOV2 extant/ancestor pair has entered preliminary laboratory work. The
YtvA and EL222 pairs remain documented here as prepared designs. These status
labels do not imply successful comparative characterization or biochemical
validation.

## Experimental Status

Experimental work on the AsLOV2 ancestral/extant comparison is ongoing and
should be considered preliminary. Expression, purification, and spectroscopy
are in progress, with repeated measurements and optimization continuing.

These activities are treated only as evidence of experimental progress. They
do not yet establish validated conclusions about comparative expression,
purification, spectroscopy, photochemical behavior, or dark-recovery kinetics.

Raw experimental files have not been uploaded, including:

- UV-visible spectra and dark-recovery measurements
- SDS-PAGE or other gel images
- purification chromatograms
- cloning, expression, and purification records
- laboratory notebook entries

Accordingly, this repository does not currently report validated experimental
results or conclusions about spectral properties, kinetics, stability,
oligomeric state, or signaling activity.

## Current Next Steps

- continue expression, purification, and spectroscopy optimization
- repeat measurements before comparative interpretation
- preserve raw spectra, gel images, chromatograms, and experimental metadata
  before adding experimental conclusions to the repository
- extend experimental work to additional selected constructs only after their
  preparation status is confirmed

## Repository Evidence

- `trees/` and `trees/comparison/`: phylogenies and topology comparisons
- `constructs/`: candidate sequences and convergence summaries
- `constructs/construct_design/`: reference-guided extraction records
- `results/construct_design/`: reviewed design rationale and final FASTAs
- `scripts/`: compact analysis utility retained for reproducibility

Large intermediate alignments, ASR state files, indel-aware working outputs,
and raw AlphaFold predictions are maintained outside this compact GitHub
repository.

## Documentation Priorities

Before publication or a formal project release:

1. Record sequence provenance and complete filtering counts.
2. Record exact software versions and commands.
3. Create a canonical candidate manifest linking node IDs to source trees.
4. Add checksums for canonical inputs and final construct FASTAs.
5. Add experimental data only with raw files, metadata, and analysis methods.
