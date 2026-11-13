# Project Status

Last documentation review: November 13, 2026

## Progress update — November 13, 2026

This update records progress reported by the researcher and the supplied October
30 analysis report. The underlying trajectories and analysis outputs were not
independently reanalyzed for this documentation update.

### MD comparison and mutant screening

The analysis now distinguishes the historical WT144 model, the rebuilt WT157
model and Anc157. Construct length, residue correspondence, starting coordinates
and simulation protocol must be recorded for each comparison.

A nine-model mutant MD set was planned. Four single-mutant systems had completed
production and analysis as of the October 30 report, each with one 10 ns
trajectory and 1,000 stored frames. Completion of the other five systems has not
been verified in this update. Model identities and numerical results remain
external under the existing public-release scope.

The four reported mutant trajectories used engine HMR with a 4 fs timestep.
Their existing WT157 reference used non-HMR dynamics with a 2.5 fs timestep.
These comparisons provide descriptive observations with a protocol confound.
Independent simulation repeats and matched-protocol comparisons remain pending.
Time blocks and stored frames are not independent simulation replicates.

### Analysis scope and experimental validation

Experimental mutation validation is in preparation. Structural metrics and
kinetic-model-derived activation quantities require further validation of their
relationship. Specific mechanistic predictions and mutation designs are retained
in internal research records.

Contact occupancy, hydrogen-bond geometry and RMSF quantify structural features
within the sampled trajectories. Their energetic contributions to activation
enthalpy or entropy have not been quantified. The fixed light-state Cys–flavin
adduct model does not sample the chemical recovery reaction or its transition
state.

The earlier sections below retain the historical project record. Numerical
comparative results, new sequences and research figures remain external according
to [the public-release scope](public_release_scope.md).



## HMR documentation update

The current HMR comparison is now documented in [hmr_analysis.md](hmr_analysis.md):
one confirmed production run per system. gREST convergence was not confirmed
for the evaluated length and metrics; see [gREST review](grest_analysis.md). Existing analysis and FPS/ChimeraX code snapshots are
preserved under `analysis/hmr/`. No analysis was rerun. Unpublished comparative
results and research figures remain withheld pending release clearance; see
[public_release_scope.md](public_release_scope.md). The older experimental
status below describes what is released here, not an assertion that no later
local measurements exist.

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

## Historical Next Steps (October 21, 2026)

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
