# Project documentation

## Purpose

Connect completed construct preparation, ongoing experimental characterization,
external MD preparation, current HMR analysis and the separate exploratory gREST
branch. This directory records methods and status; it does not contain raw
experimental data or production simulation inputs.

## Files and their roles

| Document | Role |
|---|---|
| [Project status](project_status.md) | Completed design work, experimental progress and remaining tasks |
| [Construct design](../results/construct_design/README.md) | Existing boundary rationale and final construct records |
| [HMR guide](../analysis/hmr/README.md), [HMR methods](hmr_analysis.md) | Structural-analysis entry point and definitions |
| [gREST guide](../analysis/grest/README.md), [gREST review](grest_analysis.md) | Exploratory execution and convergence evidence |
| [Software environment](software_environment.md) | Recorded versions and unresolved provenance |
| [Public release scope](public_release_scope.md) | Included code/diagnostics and excluded research materials |

## Required inputs and dependencies

Experimental interpretation requires raw spectra, recovery traces, conditions,
replicate metadata and fitting provenance outside this repository. MD
preparation requires reviewed construct-to-MD mapping, topology/coordinates,
force-field/adduct and HMR settings, equilibration/restart provenance, production
inputs and logs. Existing records show simulation preparation and production
occurred; this clone does not provide a complete preparation recipe or certify
all original inputs/builds. Package requirements depend on the chosen branch.

## Execution order / example commands

Read the [main workflow](../README.md#project-workflow-and-status), then the
relevant analysis README. There is no executable experiment or MD-preparation
pipeline in this directory. No command is supplied as a replacement for missing
private inputs. Follow the conditional commands in the HMR, frame-selection,
ChimeraX or gREST guide only after restoring the documented external workspace.

## Outputs

Documentation of methods, release scope, status and interpretation limits.
Research measurements, new numerical analyses and rendered structures are not
produced by this documentation update.

## Validation status and limitations

English prose, relative links and repository navigation are checked. Historical
validation records are distinguished from fresh checks. Experimental comparative
claims remain unreleased, HMR uses one confirmed run per system, and gREST
convergence remains unconfirmed. Contact/RMSF observations do not establish dark
recovery causality; A447 loss is not a denaturation rate or Tm, and Eyring
activation quantities are not direct evidence of structural entropy.
