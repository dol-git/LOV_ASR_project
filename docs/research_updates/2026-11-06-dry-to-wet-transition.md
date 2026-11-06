# AsLOV2 ASR — dry-to-wet handoff (public status)

**Date:** 2026-11-06  
**Phase:** exploratory computation and structural interpretation → proposed wet-lab validation  
**Release status:** methods and project-level status only. Unpublished measurements, residue-specific effect sizes, frames, figures and private coordinates remain external.

## Research question

Can evolutionary substitutions in ancestral/extant AsLOV2 change the thermal dark recovery of the covalent FMN–Cys light adduct, and can prospective mutations help distinguish candidate structural coupling from local chemistry?

## Completed research activities

- Prepared ancestral/extant AsLOV2 candidates through the documented alignment, phylogenetic, ASR and construct-design workflows.
- Pursued expression, purification and absorption-based photocycle characterization; comparative numerical results are **not released** here.
- Completed exploratory short conventional/HMR fixed-adduct simulations, contact/strict-hydrogen-bond inspection, local pocket rotamer analysis, and selected reciprocal *in silico* perturbations. Original trajectories, parameter files and output tables remain external.
- Examined candidate core–Jα interaction rearrangements using predicted-structure visualization. This is preliminary visual inspection, **not** an independent validated structural ensemble or final quantitative PDB hydrogen-bond analysis.
- Reviewed existing Eyring/Arrhenius analyses; the finalized common fitting dataset, weighting, confidence intervals and parameter covariance still require review.

## Interpretation boundaries

**Contact is not hydrogen bonding.** Historical strict HB analysis used donor–acceptor distance ≤3.0 Å and D–H–A angle ≥150°, whereas heavy-atom contact used a separate <4.5 Å definition. Spatial path existence is not a proton-transfer relay or causal allosteric mechanism.

Observed contact/rotamer Shannon entropy is bin-dependent structural-state diversity, not the absolute thermodynamic entropy of the light-adduct reactant, and not the Eyring activation entropy. State populations and dynamic coupling are not established to be converged.

Baseline systems differ in FMN torsion supplementation, hydrogen-mass repartitioning/time step and some engine provenance; all primary structural comparisons involve short single trajectories. Fixed-adduct classical MD cannot calculate adduct cleavage or its chemical transition state.

Existing adduct geometry checks do not establish a WT-specific mechanical-strain mechanism; preliminary state-population differences remain candidate explanations, not proof.

## Handoff and next verification gates

1. **Geometry audit:** compare the exact predicted-structure atom distances, side-chain orientation and HB geometry with prior trajectory definitions; verify numbering by explicit construct-to-model sequence correspondence.
2. **Sequence-source audit:** establish the actual expression plasmids, corrected backbone, tag/linker and coding junctions before designing mutant DNA or primers. Previously archived computational FASTA records alone do not establish an order-ready insert.
3. **Reciprocal substitution assay:** prepare parent and matched candidate mutants, verify their sequences and sample purity, then measure the light-adduct dark recovery in matched buffer/illumination/temperature conditions.
4. **Kinetics and thermodynamics:** preserve traces, replicate identities, raw rate fits and temperature calibration. Compare the same justified Arrhenius/Eyring temperature window and report uncertainty and parameter covariance.
5. **Interpretation:** judge whether mutation changes measured chemical recovery. Additional structural and pocket-sensitive experiments may follow, but mechanistic causality requires independent evidence.

## Provenance and public-release note

This document is a **research status update**, not an independent rerun of MD, new analysis result, structure verification, or public release of unpublished values. See [project status](../project_status.md), [HMR methods](../hmr_analysis.md), [gREST review](../grest_analysis.md) and [public-release scope](../public_release_scope.md).
