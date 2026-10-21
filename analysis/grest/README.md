# Exploratory gREST sampling review

## Purpose

Document an exploratory MD branch separate from [HMR structural analysis](../hmr/README.md).
WT and Anc each have one reviewed 50 ns exchange set of four replicas, plus an
earlier 10 ns pilot whose independence is unresolved. Convergence was not
confirmed for the evaluated length and metrics, so gREST was not used as primary
evidence for quantitative ensemble comparisons.

## Files and their roles

| File | Role |
|---|---|
| [analyze.py](Heo/analysis/grest50/analyze.py) | Inventory, DCD/log/REM QC, exchange diagnostics and physical-state conversion |
| [analyze_convergence.py](Heo/analysis/grest50/contact_distances/convergence_10ns_50ns/analyze_convergence.py) | Pilot versus 50 ns and within-state block distance comparisons |
| [Block-population analysis](Heo/analysis/grest50/clustering/block_population_stability/analyze.py) | Frozen-label fractions, JSD/TV and descriptive heterogeneity |
| [submit_50ns.example.sh](submit_50ns.example.sh) | Adapted launch recipe without personal account/path settings |
| [exchange_diagnostics.csv](exchange_diagnostics.csv) | Six existing exchange-acceptance rows, not newly calculated results |
| [code_provenance.csv](code_provenance.csv), [external_evidence_index.csv](external_evidence_index.csv) | Source/repository hashes and external evidence references |
| [Detailed review](../../docs/grest_analysis.md) | Conditions, completed lengths, state definitions and convergence limitations |

## Required inputs and dependencies

Python, NumPy, SciPy and Matplotlib are needed; conversion additionally requires
external GENESIS `remd_convert`. Exact production build/environment provenance
is incomplete. The code expects `Fugaku_Heo_20260928/`, `Heo/` and `tools/` under
an external workspace root with the original directory layout. Supply private
DCD/REM/log/restart, input/stdout, PSF/PDB/force-field files, pilot data, converted
state-1 data, frame mappings, distances and existing cluster memberships.
None of these private materials is bundled. The launch recipe also needs the
site's MPI/resource allocation and reviewed production input.

## Execution order / example commands

These are historical interfaces, only for a separate restored workspace, from
its `Heo/analysis/grest50/` directory. They are not runnable from a fresh clone:

```bash
python analyze.py prepare
python analyze.py convert WT
python analyze.py convert ANC
python analyze.py immutability
```

The distance-comparison script subsequently needs existing 50 ns distance tables
and reprocesses pilot coordinates. The block-population script needs existing
global cluster labels/populations; it does not recluster. Run each from its own
restored directory only after prerequisites are supplied. Both execute at module
top level, including on import. They can generate figures and reports; they were
not executed for this update. The launch example is an optional production
recipe, not part of the completed documentation checks and not a new job submission.

## Outputs

Historical outputs include inventory, exchange acceptance/occupancy/round-trip
records, state traces, conversion inputs, state-1 PDB/DCD and frame mappings,
coordinate QC, block distance/population summaries and figures. Private outputs
remain external. The block report keeps the existing filename `REPORT.ko.md`
while its generated prose is now English; identifiers and numerical calculations
are unchanged.

## Validation status and limitations

Two Python snapshots retain their original hashes. One has English report-prose
edits recorded with separate source/repository hashes; the numerical prefix is
unchanged. This update checks syntax, links and provenance without rerunning
analysis or conversion. Historical checks record REM/log and physical-state
mapping consistency. Parameter-space exchange mixing does not establish
structural convergence. Comparisons of different effective-temperature states
alone do not show nonconvergence. Replicas are not independent runs, the pilot
is not an established independent repeat, and Anc's tempering region differs
from its corresponding Jα analysis region. Same-state block changes and missing
independent validation limit quantitative ensemble claims.
