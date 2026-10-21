# HMR analysis: existing implementation and execution order

Reviewed 2026-10-21 by reading existing scripts and saved provenance/validation
records. No trajectory analysis, frame selection, rendering or numerical result
generation was performed for this repository update.

## Verified scope summary

This is a methods/scope table, not a comparative results table. Values were
cross-checked against source_selection.md, METHODS.md, source manifests and
script assertions in the local historical workspace.

| Item | Confirmed scope |
|---|---|
| Independent production runs used | WT run01: 1; Anc run01: 1 |
| Stored frames per system | 1,000, 0.01–10.00 ns, 10 ps interval |
| Conditions recorded | HMR, 298.15 K NPT, 4 fs, fixed Cys–FMN adduct |
| Alignment/contact core | WT11–110 ↔ Anc10–109, 100 Cα pairs |
| Operational Jα | WT120–142 ↔ Anc119–141, 23 Cα pairs |
| Contact grid | 100 × 23 = 2,300 residue pairs |
| Main contact definition | Minimum protein-heavy-atom distance <4.5 Å |
| Sensitivity windows | Full; first 500; last 500 frames, not independent runs |
| Diversity display | Global medoid + 7 farthest-point additions per system |

The current comparison excludes gREST because of convergence problems. No
separate conventional production run or additional independent HMR replicate
was confirmed. Backup copies, tutorial copies and equilibration do not increase
these counts. Production is included in full; no new equilibration cutoff is
implied.

## Alignment and RMSF

`Ja_literature_review_20261019/analyze.py` reuses the direct MD residue
correspondence and per-frame transforms from `Ja_core_visual_20261019`.
The reference is centered WT production frame0 core Cα; the transform is
`(raw - center) @ rotation`. Jα is excluded from fitting. There is no new
Kabsch fit or unwrap in this script. Prior coordinate-continuity checks are
reused. Mapping does not establish the original ASR alignment/construct source;
WT119–124 remains unresolved.

RMSF is `sqrt(mean(sum((r - mean(r))**2)))` for each matched Cα, with a
separate mean for each full/half window and the same core transforms. It is
unweighted. Historical GENESIS backbone/average-reference RMSF is a different
definition and must not be substituted for this result. Core-aligned Jα RMSD
to its own first frame and Jα self-fitted RMSD measure different quantities.

## Contact distance and controls

Atom identity is checked against PSF index/name/resname and a prior validated
identity table. Hydrogen exclusion uses element/CHARMM identity, **not HMR
mass**. Protein cysteine atoms in FMC are included; flavin, water and ions are
excluded. `distance_array` receives Å coordinates and **no box argument**;
there is no minimum-image correction at that call. Its validity relies on the
previously checked intact protein coordinates.

Minimum distances are calculated for all 2,300 pairs and all frames. The
reporting union contains pairs with distance <5 Å in at least one frame of
either system. Contacts use strict `<`, not `<=`; sensitivity cutoffs are
4.0/4.5/5.0 Å. Contact fraction is observed contact frames divided by frames,
not equilibrium population, energy, a hydrogen bond or a salt bridge.

Backbone N/CA/C/O distances are a separate geometric control.
`matched_atom_control.py` intersects heavy atoms by name, element and CHARMM
type at corresponding residues. It can remove mutation-specific terminal
functional groups; a reduced difference is not proof of absent interactions.

The descriptive screening rule requires absolute full 4.5 Å fraction difference
≥0.20, same sign and absolute difference ≥0.10 at every cutoff, and same sign
and ≥0.10 in both halves at 4.5 Å. Display pairs are ranked symmetrically by
absolute full 4.5 Å difference, ties by WT core then Jα residue number. These
rules do not test independent-run reproducibility or establish causality.

DSSP uses MDTraj `compute_dssp(simplified=False)` on an analysis-only backbone
topology, Å converted to nm. FMC cysteine is represented as CYS and histidine
variants as HIS only in that topology. Only H counts as α, not G/I. Original
PSF and coordinates remain unchanged. No p-values, SEM or bootstrap uncertainty
from treating consecutive frames as independent are justified here.

## Frame selection and ChimeraX

`Ja_core_visual_fps8_20261019/select_frames.py` reads the existing aligned
Jα cache and source transforms. It checks every raw protein frame for finite
coordinates, covalent lengths strictly between 0.5 and 3 Å, and cache agreement
<1e-5 Å. Failure stops execution rather than silently discarding frames; this is
an operational QC rule, not a physical state definition.

Distance is Jα Cα RMSD in the common-core frame without Jα refitting. The global
medoid minimizes the sum of distances. Seven subsequent selections maximize
distance to the nearest chosen frame; ties take the earliest source frame.
Selection is independent within each system using the same algorithm. Eight
frames are not eight physical states, observation frequencies or equilibrium
populations. PDB output uses actual PROA coordinates with only the existing
rigid transform, not averaging or interpolation.

`build_viewers.py` reads those PDBs and creates manifest-based camera config and
`chimerax_view16.py`. It pools both systems to set shared full/Jα orthographic
cameras. The viewer checks hashes, chain X and 23 Jα residues; commands are
`lovview WT full`, `lovview Anc full`, `lovview WT ja`, `lovview Anc ja`.
It displays the medoid core in gray, WT/Anc Jα blue/orange, and additional Jα
frames with 40% transparency. Display core (WT1–119/Anc1–118) differs from fitting
core. No MatchMaker, per-system zoom or coordinate editing is performed.
Saved native ChimeraX 1.11.1 headless checks passed; GUI rendering and GUI labels
were not verified. These are historical checks, not checks rerun here.

## Historical execution order and external inputs

The committed scripts are byte-identical historical code snapshots, not a
standalone runnable release. See `analysis/hmr/code_provenance.csv` for SHA256.
They execute work at module top level: **do not import them for inspection**.
They write into their own directories and may overwrite outputs. For any future
execution, copy them to a new, separate workspace with the original sibling
layout and input manifests; do not run in the preserved source workspace.

1. Prerequisites: validated topology/trajectory, atom identity table, direct MD
   mapping, transforms, source hashes and aligned coordinate caches. Recreate
   the `Ja_core_visual_20261019`, `Ja_sheet_tilt_20261019`,
   `Heo/analysis/reanalysis_20261006/00_provenance/atom_identity` sibling paths.
   Manifests contain external input paths; they are not included in this public
   repository. `analyze.py` also checks prior mapping and frame-selection notes.
2. `analyze.py` → `matched_atom_control.py`: distance arrays, RMSF, DSSP, contact
   fractions and matched-atom sensitivity. These scripts are included.
3. Historical downstream sequence: `derive.py` → `figures.py` → `report.py` →
   `verify.py`. Those local scripts/templates and unpublished outputs are not
   bundled; this is an execution record, not a command runnable from this clone.
   Historical figures used the earlier 5+5 frame display, not the newer FPS set.
4. FPS branch: `select_frames.py` → `build_viewers.py` → open
   `chimerax_view16.py` in ChimeraX. Also requires prior
   `Ja_core_visual_8_20261019` for protected-file checks and an existing empty
   output `final_pdb/` directory. It uses the analysis source_manifest hashes.
   Historical native viewer/algorithm validation scripts are not bundled.

The older core100 transform-generation pipeline and private inputs are not
included. Fresh-clone execution is therefore intentionally incomplete. Package
installation, MD execution and regeneration were outside this update.
