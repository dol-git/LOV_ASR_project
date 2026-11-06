# Planned reciprocal-mutant validation — AsLOV2 ASR

**Date:** 2026-11-06  
**Status:** proposed; sequence/source approval and ordering not completed.  
**Public scope:** mutation-site identities and unpublished structure/kinetics effect sizes are held in the controlled lab record pending release clearance.

## Question and design

Test whether perturbing a candidate core–Jα interface interaction affects dark recovery of the FMN–Cys light adduct. Compare **WT and Anc parents** with **one matched reciprocal substitution in each corresponding background**. A shift in recovery rate would motivate—but would not by itself prove—a particular structural or chemical mechanism.

The site is tracked against the actual expression construct and correspondence mapping in the lab-held record. No external sequence, DNA FASTA, primer or ordering specification is certified by this GitHub document.

## Sequence and ordering gates

- [ ] Confirm the finalized WT and ancestor expression-plasmid DNA (including any corrected backbone).
- [ ] Reconcile the construct, predicted-structure and MD residue mappings by sequence alignment rather than a fixed offset.
- [ ] Validate both parent and variant translations, codons, tags, linkers, cloning junctions and intended single substitutions.
- [ ] Check that no unintended amino acid or coding changes are present, and archive source hashes and a sequence diff.
- [ ] Decide between site-directed mutagenesis and gene synthesis only after verified plasmid and insert information is available.
- [ ] Mark constructs **NOT ORDER READY** until these checks pass and the supervisor approves the sequence.

## Matched experimental readouts

- Confirm expression, FMN incorporation, purification quality and consistent sample preparation across parent and variant proteins.
- Remove residual imidazole; match buffer composition, pH, ionic strength and illumination protocol.
- Record dark/illuminated UV–Vis signatures, time-dependent recovery traces, actual sample temperature, fit residuals and independent replicate IDs.
- Analyze recovery without forcing a one-exponential model if the raw traces disagree.
- Choose a shared temperature window appropriate to the kinetic model; scrutinize high-temperature deviation separately.
- Report fitted k and tau, Arrhenius/Eyring activation quantities, their uncertainties/covariance and temperature-window sensitivity.

## Distinguishing outcomes

| Outcome | Limited inference |
| --- | --- |
| Reproducible kinetic change with mutation | The candidate site contributes to measured recovery; molecular mechanism remains open. |
| Local structural change but no kinetic change | Weakens a simple contact-controls-chemical-rate model. |
| Background-dependent mutant effect | Motivates further epistasis/structural tests, after checking construct and measurement differences. |
| No detectable mutant effects | Weakens this specific structural causal explanation without excluding other pocket routes. |

## Boundaries

Predicted structures and short, non-identically parameterized fixed-adduct trajectories are exploratory evidence. Contact occupancy, local Shannon entropy and adduct geometry cannot be equated to activation entropy, reaction barrier or transmitted mechanical tension.

No sequence has been released as order-ready, no variant has been reported experimentally validated, and no synthesis order has been submitted under this plan.

See [handoff status](../research_updates/2026-11-06-dry-to-wet-transition.md) and [release policy](../public_release_scope.md).
