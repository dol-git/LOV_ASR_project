# Repository Organization

Reviewed November 13, 2026. The review inspected all 127 previously tracked text
files, including code, sequence/alignment records, tables and guides.
The [file inventory](repository_inventory.csv) records retention and removal
categories, including newly added guides.

| Category | Decision |
|---|---|
| Methods, source code, requirements and provenance | Preserve stage-specific definitions and original calculations |
| Scientific sequences, alignments and comparison tables | Retain existing records; distinguish final constructs from intermediate designs |
| Identical tree duplicates | Remove eight duplicate files; canonical files remain directly under trees/ |
| Empty FASTA placeholders | Remove four empty files; populated selection records remain |
| Stage-specific duplicate inputs | Retain when expected by scripts or useful for provenance |
| Historical node metadata | Retain with an explicit source-mapping caveat |
| Main README | Condense to project scope, status, navigation and construct overview |
| Progress and release documents | Consolidate status and content policy; remove procurement and editing-session narrative |

## Records worth adding when verified

- Source accession lists, filtering counts and database releases.
- A construct manifest linking sequence checksums, source trees, node IDs and boundaries.
- Run manifests recording construct, starting structure, protocol, seed, force field and software build.
- A clear external-input index for each analysis snapshot.

Add verified records as they become available. Existing historical scripts depend
on external inputs and are preserved at their original paths. This organization
review does not establish release clearance for pre-existing sequences or rerun
scientific analyses. Code execution would require the external research workspace.
