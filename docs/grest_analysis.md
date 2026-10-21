# gREST execution and convergence review

Reviewed 2026-10-21 against existing inputs, production logs, analysis code,
reports and diagnostic CSVs. No analysis was rerun and no DCD was reprocessed.
Completed lengths were checked against the final step in every replica log,
not just requested input lengths. Structural evidence was reviewed from saved
reports; this does not constitute an independent reproduction.

**Convergence was not confirmed for the evaluated simulation length and metrics,
so these calculations were not used as primary evidence for quantitative
ensemble comparisons.**

## Execution scope and conditions

| Item | WT | Anc |
|---|---|---|
| Reviewed 50 ns production sets | 1 | 1 |
| Exchanging execution replicas per set | 4 | 4 |
| Final completed step per replica | 12,500,000 | 12,500,000 |
| Timestep / nominal completed length | 0.004 ps / 50 ns | Same |
| Saved frames per replica | 5,000, 0.01–50.00 ns | Same |
| Bath / thermostat | NVT, 298.15 K, BUSSI | Same |
| Solute parameter states 1–4 | 298.15 / 331 / 365 / 400 K | Same |
| REST settings | 1 dimension, `param_type1=D L` | Same |
| Actual production tempering selection | PROA residues 120–142 | PROA residues 120–142 |
| Corresponding analysis Jα region | 120–142 | 119–141 |

The inputs specify CHARMM, PME, vdW force switching, switch/cutoff/pairlist
radii of 10/12/13.5 Å, VVER, rigid bonds and PBC. HMR settings are
`hydrogen_mr=YES`, `hmr_ratio=3.0`, `hmr_ratio_xh1=2.0` and
`hmr_target=solute`. These are recorded input settings, not verification of all
force-field files, the production binary or physical accuracy. Coordinates and
exchange attempts are scheduled every 2,500 steps (10 ps). All eight 50 ns logs
end at INFO step 12,500,000 with printed time 50,000.0024 ps. The nominal time
axis uses step × 0.004 ps; rounding differences are not additional sampling.

An earlier 10 ns pilot also has four replicas per system. All eight pilot logs
reach 2,500,000 steps (printed time 10,000.0005 ps); existing manifests and QC
record 1,000 frames per replica at 10 ps intervals. The pilot distance traces
are not identical to the first 10 ns of the 50 ns set. Independence of restart
histories and initial conditions has not been established: the pilot is neither
counted as an additional independent repeat nor concatenated into a 60 ns run.
Its complete production settings are not assumed identical to the 50 ns set.
Short tests, equilibration and archive copies do not add independent runs.

## Replicas, walkers and the physical state

- **Execution replica / walker:** follow a raw replica DCD together with its REM
  history. In this recorded implementation, coordinates are written before
  parameter assignments are updated. A walker's Hamiltonian state changes with
  exchanges. Four interacting walkers form one exchange set, not four
  independent production runs.
- **Parameter state:** a fixed REST condition. Higher K labels describe the
  effective solute tempering parameters for selected terms, not separate NVT
  simulations with the entire bath heated to those temperatures. `D L` tempers
  selected dihedral and Lennard–Jones terms.
- **Physical state:** `remd_convert` with `convert_type=PARAMETER` and
  `convert_ids=1` selects the replica holding state 1 at each saved time. State 1
  is the 298.15 K baseline. GENESIS source defines
  `coeff_full=T_bath/T_rest` and `coeff_half=sqrt(coeff_full)`, giving unity
  scaling at state 1. This follows from source and input settings, not a
  separately printed runtime scaling check.

Saved frame i (1-based), at step i×2500, uses the state in REM step
(i−1)×2500. Existing QC checked REM/log transition consistency and exact
agreement of converted coordinates and boxes with the source state-1 replica
for every frame. The later walker audit also checked state-1 owner mapping.
The converted DCD's default time header does not represent the true time axis;
use the 10 ps/frame recorded in `physical_frame_mapping.csv`.

A physical-state trajectory joins coordinates from different walkers into
parameter-sorted samples. It does not represent uninterrupted single-molecule
kinetics, transition rates or lifetimes. Matching a high-state geometry to a
later state-1 geometry of the same walker is a provenance/geometry diagnostic,
not a physical return time or evidence of structural equilibrium.

## Confirmed sampling diagnostics

[exchange_diagnostics.csv](../analysis/grest/exchange_diagnostics.csv) contains
six rows copied from existing `qc/{WT,ANC}/exchange_acceptance.csv` files.
Acceptance rates and exchange events were not recalculated for this review.

| System | State 1↔2 | State 2↔3 | State 3↔4 | Attempts per pair |
|---|---:|---:|---:|---:|
| WT | 24.68% | 25.84% | 28.16% | 2,500 |
| Anc | 21.68% | 25.36% | 28.56% | 2,500 |

The original code deduplicates both replicas' exchange records by
`(step, state pair)` and checks their agreement. Existing summaries record zero
REM/log mismatches and 212 WT / 211 Anc completed, non-overlapping A→B→A trips
between the endpoint states. Walker occupancy tables and parameter traces also
record movement along the ladder. Incomplete final trips are censored
separately; the final post-exchange point adds no residence time. These checks
show **parameter-space mixing**, not sufficient sampling of slow structural
transitions or equilibrium. Structural effect sizes, sequences and coordinates
are not released here.

## Structural metrics and comparison targets

The primary convergence question is whether distributions reproduce across
time windows **within the same system at physical state 1**. Distribution
differences between effective-temperature states alone are not evidence of
nonconvergence: their target Hamiltonians differ. Pooling all states without
appropriate reweighting does not yield the baseline ensemble.

| Existing check | Definition and comparison | Supported interpretation |
|---|---|---|
| Jα RMSD | Unweighted core-backbone fitting; Jα N/CA/C/O relative to its own reference, without Jα refitting | Block variation and fitting sensitivity |
| RMSF | Full-window and block fluctuation/profile comparisons | Window dependence; a single average cannot establish convergence |
| Contact distances | Minimum heavy-atom distances for seven corresponding residue pairs; median/IQR/ECDF in five 10 ns blocks | Some pairs remain block dependent, not a universal failure of all pairs |
| Pilot versus 50 ns | Distance distributions with matching pair/atom definitions | Comparison of observed datasets, not an independent-sample test |
| Fixed cluster fractions | Global state-1 labels retained across blocks at 1.00/1.25 Å cutoffs | Descriptive population stability |
| Extended descriptor audit | Contact counts, COM, Jα orientation, DSSP H fraction and PCA axis quality across blocks/cumulative samples | Descriptor-specific drift or block sensitivity; a stable label does not establish equilibrium |

The 50 ns distance report uses five blocks, (0,10] through (40,50] ns.
ECDF difference is sup|F−G|; JSD is base-2 divergence in bits, not its square-root
distance; TV is ½Σ|p−q|. WT/Anc cluster labels are defined separately, so no
cross-system cluster-vector JSD is calculated. The full distribution contains
each compared block, and cumulative samples overlap. Apparent improvement
against the full distribution alone is insufficient; block comparisons take
priority. No p-values/SEM assuming independent consecutive frames or
predeclared convergence pass threshold were provided.

Existing reports retain shifts in the location/width of some state-1 distance
distributions and block dependence of fixed-label fractions. Some contact and
structural metrics reproduce more closely, so these findings do not establish
nonconvergence of every observable. Together with the lack of independent
repeats, they support the conservative decision not to use these samples as
primary evidence for stable quantitative ensemble comparisons. They are not a
formal statistical proof of nonconvergence based merely on differing distributions.

## Assessment of earlier conclusions and unresolved issues

1. The earlier phrase “convergence problems” lacked scope. The conclusion is now
   restricted to convergence not being confirmed for the evaluated length and
   metrics. This does not establish convergence of the HMR analysis either.
2. Initial fitting cores were WT1–119/Anc1–118. After checking N-terminal
   sensitivity, stable-core fitting used WT11–119/Anc11–118; homologous-core
   sensitivity used WT12–119/Anc11–118. These differ from the recent HMR core100
   Cα definition. An initial core-RMSD high tail alone cannot establish global
   core collapse or denaturation; statistics with different definitions must not
   be combined.
3. Anc production tempering residues 120–142 differ from the corresponding
   Jα119–141 analysis region, as confirmed in input and stdout. Changing a
   postprocessing selection does not correct the production Hamiltonian. This
   is a separate comparison-condition limitation; its effect size is unresolved.
4. The pilot README's default Jα120–142 description differs from later
   corresponding-region analysis at 119–141. The initial 50 ns RMSF report lacks
   DSSP, whereas a later separate reanalysis includes MDTraj DSSP. Reports from
   different stages do not describe one identical analysis version.
5. Effective independent sample size, uncertainty, independent-run
   reproducibility and rare-state stability remain unverified. Exchange rates,
   round trips, same-state walker differences or owner fractions alone do not
   establish success or failure. The fixed adduct does not sample dark recovery
   or its transition state; cluster Shannon H is not thermodynamic structural
   entropy or direct evidence of Eyring ΔS‡.
6. Future work requires corresponding tempering-region and input/force-field
   provenance checks, exchange sets with independent initial conditions, and
   predefined same-state metrics with correlation-aware uncertainty. No new
   simulation, extension, reweighting, FES or analysis run was performed here.

## Code, external materials and execution order

Start with the [gREST directory guide](../analysis/grest/README.md).
[Code provenance](../analysis/grest/code_provenance.csv) records original source
SHA256 values. Two Python snapshots remain byte-identical; the block-population
script has only report-language edits, with a separate repository SHA256 and
adaptation note. No analysis code was executed during this documentation update.
The [launch example](../analysis/grest/submit_50ns.example.sh) removes personal
paths, allocation accounts and site modules from the historical WT/Anc scripts
and takes executable/input paths from environment variables. It is not an
identical copy of a historical submission script and has not been job-tested.
Private production inputs and coordinates are required; a fresh clone is not a
complete executable workflow.

Snapshots remain under `analysis/grest/Heo/analysis/grest50/`. Copy this sibling
layout to a separate external workspace with the required source material before
any future execution. The following describes code behavior, not commands run
for this update:

- `analyze.py prepare`: source inventory, DCD record/frame checks, log/REM
  exchange QC, traces, round trips and conversion-input generation; outputs may
  be overwritten.
- `analyze.py convert WT` / `convert ANC`: invoke external GENESIS
  `remd_convert` and verify converted coordinates. `analyze.py immutability`
  checks raw files against the saved inventory.
- `contact_distances/convergence_10ns_50ns/analyze_convergence.py`: uses existing
  50 ns distance CSVs plus pilot DCD/REM/PSF, recalculates pilot distances and
  creates figures. It is not merely a CSV summary command. Most outputs use
  exclusive creation.
- `clustering/block_population_stability/analyze.py`: uses existing memberships,
  populations and frame mappings to calculate fractions/JSD/TV. It does not
  recluster, but generates reports and figures. Work runs at module top level,
  so importing it also executes analysis.

External requirements include production DCD/REM/log/restart, inputs/stdout in
`Fugaku_Heo_20260928/HMR_{WT,ANC}/`, matching `input/` PSF/PDB and force-field
files, `tools/genesis-2.1.6/bin/remd_convert`, physical-state conversions and
mappings, existing distance/atom-selection tables, cluster memberships and
populations, and the pilot analysis/source files. A complete dependency bundle
or environment export is not included. Code uses NumPy, SciPy, Matplotlib and
Python's standard library; the exact production build and complete environment
remain subject to external provenance checks.

The [external evidence index](../analysis/grest/external_evidence_index.csv)
records workspace-relative paths and SHA256 values for reviewed reports and
metadata. Main evidence includes `Heo/analysis/grest50/REPORT.ko.md`,
`FINAL_EVIDENCE_AUDIT.ko.md`, `walker_provenance/REPORT.ko.md`,
`contact_distances/convergence_10ns_50ns/REPORT.ko.md`,
`clustering/block_population_stability/REPORT.ko.md` and separate
`reanalysis_20261006/04_phase3A_ensemble/` documentation. Production inputs and
completion steps in all eight logs were read directly but are not bundled.
Uncleared sequences, coordinates, production inputs, structural numerical
results, original reports, complete logs, DCD/restart/checkpoint files,
environment folders, personal paths and account information remain excluded.
Public numbers are limited to the requested execution conditions and compact
exchange sampling diagnostics. Paths and data identifiers retain their original
spelling, including `REPORT.ko.md`.
