# HMR structural analysis

## Purpose

Describe fluctuations and Jα–core geometry in the current WT/Anc HMR comparison.
There is one confirmed production run per system: 1,000 frames over
0.01–10.00 ns. This guide covers analysis; [frame selection](Ja_core_visual_fps8_20261019/README.md)
and [ChimeraX display](chimerax/README.md) are separate stages.

## Files and their roles

| File | Role |
|---|---|
| [analyze.py](Ja_literature_review_20261019/analyze.py) | Core-aligned Cα RMSF, protein-heavy/backbone contact distances and eight-state DSSP |
| [matched_atom_control.py](Ja_literature_review_20261019/matched_atom_control.py) | Contact control using shared atom name, element and CHARMM type |
| [requirements.txt](Ja_literature_review_20261019/requirements.txt) | Recorded package versions |
| [code_provenance.csv](code_provenance.csv) | SHA256 of unchanged historical snapshots |
| [Detailed methods](../../docs/hmr_analysis.md) | Definitions, windows, mapping caveats and external input layout |

## Required inputs and dependencies

Requires the original PSF/production DCD, validated atom identity table,
`Ja_core_visual_20261019` correspondence/transforms/frame manifests/PDBs and
aligned Jα caches, plus `Ja_sheet_tilt_20261019` validation and descriptor files.
The code expects the original sibling directory layout and paths stored in the
manifests. These inputs are not bundled. NumPy, MDAnalysis and MDTraj are needed;
see [version records](../../docs/software_environment.md). The script prepends a
local `software/` directory for MDTraj; vendor environments are excluded.

## Execution order / example commands

These are historical commands for a separate workspace with all required inputs
restored, run from its `Ja_literature_review_20261019/` directory. They will not
work from a fresh clone and can overwrite files beside the scripts.

```bash
python analyze.py
python matched_atom_control.py
```

Both scripts execute analysis at module top level, including on import.
Historical downstream `derive.py`, `figures.py`, `report.py` and `verify.py` are
external and not callable from this clone. Frame selection consumes source hashes
from this analysis but does not replace full-trajectory statistics.

## Outputs

`*_frame_arrays.npz`, `*_DSSP8.npy`, `common_core_aligned_RMSF.csv`,
`contact_fractions.csv`, `contact_atom_selections.csv`, `source_manifest.csv`,
`union_pair_distances.csv.gz`, validation/summary JSON and matched-atom control
tables. The scripts do not write PDBs. These research outputs are not bundled.

## Validation status and limitations

Snapshots match their recorded source hashes; syntax is checked without running
analysis. Historical saved checks validated source hashes, coordinates and
selected contact/RMSF calculations; they were not rerun here. Current contacts
use strict distance <4.5 Å, not energy or equilibrium population. Hydrogen
selection does not use HMR mass. One run per system and two trajectory halves
cannot establish run-to-run reproducibility or convergence. Mapping at the
WT119–124 boundary and original ASR provenance remain unresolved. Contact/RMSF
differences do not establish the cause of dark recovery.
