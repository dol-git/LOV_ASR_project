# Analysis directory guide

## Purpose

Navigate the existing MD analysis snapshots. HMR structural analysis and the
exploratory gREST sampling review are separate branches from MD preparation;
gREST is not an input stage required to run HMR.

## Files and their roles

| Guide | Role |
|---|---|
| [HMR](hmr/README.md) | Common-core RMSF, contact geometry and DSSP |
| [Frame selection](hmr/Ja_core_visual_fps8_20261019/README.md) | Medoid and farthest-point examples for display |
| [ChimeraX visualization](hmr/chimerax/README.md) | Manifest-based display with shared cameras |
| [gREST](grest/README.md) | Exchange QC, physical-state reconstruction and convergence review |

## Required inputs and dependencies

Private PSF/PDB/DCD inputs, mappings and provenance caches are external. See the
[project documentation](../docs/README.md) for preparation boundaries and
[software records](../docs/software_environment.md) for confirmed versions.
The snapshot collection is not a runnable release from a fresh clone.

## Execution order / example commands

Read the relevant guide before restoring any external inputs. The historical
HMR sequence is analysis → matched-atom control → frame selection → viewer
configuration → ChimeraX. gREST has its own prepare/convert/QC and structural
stability checks. Do not import analysis scripts for inspection: several execute
work at module top level. No analysis commands were run for this documentation update.

## Outputs

Historical products include coordinate arrays, distance/RMSF tables, frame
manifests, PDB examples and diagnostic reports. Only code, provenance and the
cleared compact exchange diagnostic are included here; private products remain external.

## Validation status and limitations

Documentation checks cover English prose, relative links, directory layout,
syntax and source provenance. They do not repeat historical analysis validation.
HMR has one confirmed production run per system. gREST convergence was not
confirmed for the evaluated length and metrics; exchanging replicas are not
independent runs. See [interpretation limits](../README.md#current-interpretation-and-validation-plan).
