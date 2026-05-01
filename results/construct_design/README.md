# Construct Design for LOV Ancestral Expression Candidates

This directory contains clean construct-design records for selected ancestral LOV-domain candidates and their corresponding wild-type reference constructs.

The working files used for reference-based extraction are stored separately in `constructs/construct_design/`.

## Purpose

The purpose of this step was to define experimentally testable protein constructs that contain not only the conserved LOV core but also the architecture-specific C-terminal helix or linker region relevant to signal transmission.

## Architecture classes

Three LOV-domain architecture classes were considered:

1. Plant-like LOV2 / AsLOV2
2. LOV-STAS / YtvA
3. LOV-HTH / EL222

For each class, a wild-type reference construct and an ancestral construct were prepared.

## Final construct set

| Construct name | Reference boundary | Length | Description |
|---|---:|---:|---|
| `WT_AsLOV2_404_560` | AsLOV2 404-560 | 157 aa | Wild-type LOV2 core plus Jalpha-containing C-terminal region |
| `AncPlant_Node62_AsLOV2_404_560eq` | AsLOV2 404-560 equivalent | 157 aa | Plant-like ancestral LOV2 construct with Jalpha-containing C-terminal region |
| `WT_YtvA_20_147` | YtvA 20-147 | 128 aa | Wild-type LOV core plus LOV-STAS linker / Jalpha-like segment |
| `AncSTAS_Node227_YtvA_20_147eq` | YtvA 20-147 equivalent | 128 aa | Ancestral LOV-STAS construct with LOV-STAS linker / Jalpha-like segment |
| `WT_EL222_1_163` | EL222 1-163 | 163 aa | Wild-type LOV domain plus connector Jalpha-helix |
| `AncHTH_Node83_EL222_1_170eq` | EL222 1-170 equivalent | 148 aa | Ancestral LOV-HTH construct with Jalpha-like connector segment |

## Design rationale

### Plant-like LOV2 / AsLOV2

The AsLOV2 construct was defined using residues 404-560.

This region contains the LOV2 core and the C-terminal Jalpha-containing region. The ancestral plant-like construct was extracted from `plant_linsi_Node62` using the same AsLOV2 404-560 reference coordinates.

The ancestral sequence retained the C-terminal segment corresponding to the AsLOV2 Jalpha-containing region, so no additional extension was applied.

### LOV-STAS / YtvA

The YtvA construct was defined using residues 20-147.

This region contains the LOV core and the downstream LOV-STAS linker / Jalpha-like segment, while excluding the STAS domain body.

Alignment of the ancestral `Node227` sequence against the YtvA reference showed that the C-terminal LOV-STAS linker region was well preserved with minimal gaps. Therefore, the YtvA 20-147-equivalent ancestral construct was retained as the final LOV-STAS candidate.

### LOV-HTH / EL222

The reported EL222 domain boundary was used as the reference:

- LOV domain: residues 1-144
- connector Jalpha-helix: residues 145-163
- HTH domain: residues 164-222

Therefore, the wild-type EL222 construct was defined as residues 1-163, including the LOV domain and the connector Jalpha-helix.

For the ancestral HTH candidate, the strict EL222 1-163-equivalent extraction retained a relatively short C-terminal connector segment because of alignment gaps. Therefore, EL222 1-163, 1-170, and 1-175 equivalent ancestral constructs were compared.

The EL222 1-170-equivalent construct was selected as the primary ancestral HTH candidate because it preserved the Jalpha-like charged/helical connector segment more completely while minimizing extension into the HTH domain.

## Important note on reference boundaries

Reference-coordinate extraction does not always produce ancestral constructs with the same length as the wild-type reference. This is because the ancestral sequence may contain gaps relative to the reference alignment.

Therefore, construct boundaries were selected using two criteria:

1. Literature-supported reference boundaries.
2. Actual retained C-terminal residues in the ancestral sequence.

This was especially important for the EL222-like HTH candidate.

## Files

Main files in this directory:

- `results/construct_design/README.md`
- `results/construct_design/final_order_fastas/all_WT_and_ancestral_order_constructs.fasta`
- `results/construct_design/final_order_fastas/order_YtvA_EL222_WT_Ancestor_constructs_renamed.fasta`

Working files used to generate these records:

- `constructs/construct_design/alignments/`
- `constructs/construct_design/inputs/`
- `constructs/construct_design/outputs/`
- `constructs/construct_design/scripts/extract_by_ref_range.py`

## Current status

The following constructs have already been submitted for ordering:

- `WT_AsLOV2_404_560`
- `AncPlant_Node62_AsLOV2_404_560eq`

The following constructs were prepared as the next ordering candidates:

- `WT_YtvA_20_147`
- `AncSTAS_Node227_YtvA_20_147eq`
- `WT_EL222_1_163`
- `AncHTH_Node83_EL222_1_170eq`

## Next steps

1. Submit the YtvA and EL222 WT/ancestor constructs for ordering.
2. Read papers on AsLOV2, YtvA, and EL222 signal-transmission helices/linkers.
3. Prepare structure-prediction inputs for the final construct set.
