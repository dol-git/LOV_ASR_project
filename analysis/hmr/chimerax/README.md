# ChimeraX visualization

## Purpose

Display HMR-selected frames with shared WT/Anc cameras and distinguish the
central medoid from additional examples. This is a documentation directory;
viewer code remains in its original [frame-selection directory](../Ja_core_visual_fps8_20261019/README.md).

## Files and their roles

| File in `../Ja_core_visual_fps8_20261019/` | Role |
|---|---|
| [build_viewers.py](../Ja_core_visual_fps8_20261019/build_viewers.py) | Validates PDB selections and generates shared display configuration |
| [chimerax_view16.py](../Ja_core_visual_fps8_20261019/chimerax_view16.py) | Registers `lovview`, checks hashes and displays structures |
| `selected_frames.csv`, `chimerax_display_config.json`, `final_pdb/` (external) | Required adjacent manifest, camera settings and structures |

## Required inputs and dependencies

Viewer preparation requires NumPy and MDAnalysis. Viewing requires ChimeraX
1.11.1 with its bundled Python/NumPy and native `chimerax` modules; this script
is not a normal `python chimerax_view16.py` program. Restore the manifest,
configuration and selected PDBs beside the viewer. Their public-release status
has not been established, so these inputs are external.

## Execution order / example commands

Complete [frame selection](../Ja_core_visual_fps8_20261019/README.md), then run
`python build_viewers.py` from its restored workspace. In ChimeraX, use
**File → Open** to open that workspace's `chimerax_view16.py`. The following
commands then switch views:

```text
lovview WT full
lovview Anc full
lovview WT ja
lovview Anc ja
```

Omitting the final mode preserves the current preset. Opening the committed
viewer without its private adjacent inputs will fail; do not move the code
into this documentation directory.

## Outputs

Viewer preparation creates camera/configuration and residue-audit files.
The viewer opens 16 models and registers an interactive command; it does not
save an image or alter coordinates, bonds or model transforms. Only the medoid
core is shown in gray; WT/Anc Jα is blue/orange and additional Jα examples have
40% transparency. Full/Jα presets share WT/Anc scale and camera settings.

## Validation status and limitations

Historical ChimeraX 1.11.1 headless checks covered loading, all four switches,
camera/scale/center agreement, coordinates and ribbon visibility/transparency.
GUI rendering and labels were not verified. Those checks were not rerun here.
Display core WT1–119/Anc1–118 differs from the fitting core. No MatchMaker,
per-system fitting or automatic zoom is applied. Selected frames emphasize
diversity and do not represent observation frequencies or equilibrium populations.
