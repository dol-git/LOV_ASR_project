# Public release scope — 2026-10-21

GitHub repository metadata confirmed `dol-git/LOV_ASR_project` is public.
Only methods, code snapshots without embedded sequences/results, version
records, code checksums and a checked scope table are added in this update.

Unpublished experimental values, contact/RMSF effect sizes, residue-specific
result tables, selected frame indices, structures and all existing research
figures are withheld: their public-release status could not be established.
No new figure is generated. No existing figure was found with confirmed release
clearance, so no scientific image is included in this commit. The README's
existing workflow diagram remains available. A public results figure requires
release clearance, traceable source data and an appropriate descriptive caption.

DCD, restart/checkpoint files, raw/large arrays, PSF/PDB inputs, raw spectra,
sequences not already tracked, personal paths from manifests, credentials and
environment/vendor folders are not added. Existing tracked files and history
are retained. This scoped review does not certify previously tracked sequences
as cleared for publication.

The pre-existing untracked `docs/dataset_construct_map.md` is left untouched
and excluded from staging. No push is performed.

## Subsequent gREST documentation update

The gREST update adds reviewed methods/code, an account-free launch example,
execution conditions and a six-row exchange diagnostic table specifically within
the requested sampling review. Unpublished structural effect sizes, coordinates,
sequences, production inputs and raw reports/logs remain external. See
[grest_analysis.md](grest_analysis.md). The no-push statement above describes the
original local HMR documentation update, not the subsequent authorized pushes.

## Documentation organization and language

Public documentation and generated report prose are maintained in English;
paths and data identifiers keep their original spelling. The block-population
script's report strings were translated without changing numerical logic.
Its original source SHA256 is retained alongside the current repository SHA256
and an adaptation note in `analysis/grest/code_provenance.csv`. Other snapshots
remain unchanged. Added directory guides document external-input requirements;
no analysis or private result release accompanies this update.
