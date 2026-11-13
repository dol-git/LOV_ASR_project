# Project Status

Updated November 13, 2026.

| Area | Documented progress | Remaining work |
|---|---|---|
| ASR | Four MAFFT strategies, phylogenetic comparison and indel-aware reconstruction completed | Consolidate source accessions, filtering counts and database releases |
| Construct design | Six extant/ancestral constructs prepared across AsLOV2, YtvA and EL222 | Maintain construct-to-tree and reference-boundary correspondence |
| Experiments | Expression, purification and spectroscopy pursued for the AsLOV2 pair; mutation validation in preparation | Controlled mutant characterization and replicate measurements |
| Historical HMR | One WT and one Anc 10 ns production trajectory analyzed | Independent repeats and mapping validation |
| WT157 and mutant MD | Rebuilt WT157 comparison and mutant screening initiated | Matched-protocol comparisons and repeat sampling |
| gREST | Four-replica 50 ns exchange set per system and earlier 10 ns pilot reviewed | Structural convergence assessment |

## Mutant MD record

The October 30, 2026 analysis report documents four completed systems from a
planned nine-model set. Each contains one 10 ns production trajectory and 1,000
stored frames. Later completion status has not been established in this record.
The source report was reviewed; underlying trajectories were not reanalyzed for
this documentation revision.

The reported mutant runs used engine HMR with a 4 fs timestep. Their existing
WT157 reference used non-HMR dynamics with a 2.5 fs timestep. Protocol differences
and single-run sampling limit attribution of observed differences to mutations.
Historical WT144, rebuilt WT157 and Anc157 must be distinguished by construct,
residue mapping, starting coordinates and simulation protocol.

## Interpretation and evidence

Structural metrics describe fluctuations, contacts and geometry in the sampled
trajectories. The energetic contribution of those features to activation
enthalpy or entropy remains unquantified. Fixed Cys–flavin adduct simulations
sample the light-state structure without chemical cleavage or a reaction
transition state. Time blocks and stored frames are correlated observations.

The gREST analysis evaluates exchange mixing and structural convergence
separately. Convergence was unconfirmed for the evaluated metrics and lengths.
Detailed evidence and historical execution records are in the
[HMR methods](hmr_analysis.md) and [gREST review](grest_analysis.md).

Experimental values, mutation-specific results and detailed mechanistic
predictions are maintained in internal research records under the
[public content scope](public_release_scope.md).
