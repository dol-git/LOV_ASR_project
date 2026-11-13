# Candidate and Construct Records

| Path | Role |
|---|---|
| plant/ | Plant candidates from alternative alignments and combined input |
| bacterial/ | Bacterial candidates grouped by output architecture |
| consensus/ | Stage-specific inputs, alignments and consensus outputs |
| construct_design/ | Reference sequences, alignments, extraction script and intermediate constructs |
| metadata/selected_nodes.tsv | Historical node-selection list; source provenance requires reconciliation |
| *_selected.fasta, *_targets.txt | Historical candidate-selection records |
| *_pairwise.tsv, *_pairwise.summary.tsv | Pairwise identity and summaries |
| [Final construct records](../results/construct_design/README.md) | Canonical experimental construct rationale and final FASTA set |

## Historical records and identity definitions

`metadata/selected_nodes.tsv` contains earlier node IDs that differ from the
final design set. Use the final construct guide for current selections and
interpret each historical node ID alongside its source tree. Intermediate
extractions can differ from the final expression constructs.

`pairwise_identity.py` uses a global alignment with zero gap penalties and
reports identity across paired nongap positions. `pairwise_identity_exact.py`
uses configurable gap penalties and also reports alignment coverage. Retain both
method records and check the producing command before interpreting a table.

Some stage-specific FASTA inputs have identical bytes to upstream candidates.
They are retained to preserve script input paths and provenance. Four empty
FASTA placeholders in the former `fasta/` directory were removed.
