# HMR frame selection

## Purpose

Choose actual frames that illustrate Jα structural diversity after common-core
alignment. Each system contributes a global medoid and seven farthest-point
additions. This display selection does not estimate populations or replace the
full-trajectory [HMR analysis](../README.md).

## Files and their roles

| File | Role |
|---|---|
| [select_frames.py](select_frames.py) | Source/QC checks, deterministic medoid/FPS selection and 16 final PDBs |
| [build_viewers.py](build_viewers.py) | Reads selected PDBs, computes pooled camera presets and writes viewer/config |
| [chimerax_view16.py](chimerax_view16.py) | Opens the selected structures in ChimeraX |
| [Visualization guide](../chimerax/README.md) | Display setup, commands and validation limits |

## Required inputs and dependencies

Requires NumPy, SciPy and MDAnalysis, original PROA PSF/DCD, the sibling
`Ja_core_visual_20261019/` manifests, correspondence, transforms and aligned
`*_common_core_Ja_CA.npy` caches, plus the HMR analysis `source_manifest.csv`.
Prior `Ja_core_visual_8_20261019/` files are included in protection checks.
Paths inside manifests must resolve. Restore inputs in a separate workspace;
create an empty `final_pdb/` output directory there before running. Private
coordinates, caches and manifests are not included in the clone.

## Execution order / example commands

After HMR source hashes are available, run from the restored
`Ja_core_visual_fps8_20261019/` directory:

```bash
python select_frames.py
python build_viewers.py
```

These commands are conditional on the external inputs, not ready-to-run clone
instructions. Both scripts work at module top level and write beside themselves;
use a new output workspace. The second command prepares the
[ChimeraX stage](../chimerax/README.md).

## Outputs

`final_pdb/` contains eight real frames per system; `selected_frames.csv` records
rank, source frame/time and role. Other products include
`raw_all_frame_QC.csv`, `selected_pairwise_distances.csv`, selection scores,
source/validation tables and verification JSON. Viewer preparation writes
`chimerax_display_config.json`, `chimerax_actual_residues.csv`,
`viewer_template.py.txt` and `chimerax_view16.py`. No products were regenerated here.

## Validation status and limitations

The selection code first checks all protein frames for finite coordinates,
covalent bond continuity and cache agreement, and stops on failure. It applies
existing rigid transforms without averaging or Jα refitting. Saved historical
checks cover medoid/FPS steps, tie-breaking, coordinate preservation and PDB
rounding; this update checks syntax and provenance only. Earliest source index
breaks ties. Diversity selection can emphasize rare arrangements; eight frames
are not eight physical states, frequencies or equilibrium populations. Only
HMR run01 per system is used, and gREST is excluded from this selection.
