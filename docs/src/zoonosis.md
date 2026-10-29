# Zoonosis study and migration

The framework's generic inputs are reaction systems with species/parameter maps,
`TreeProblem`, `ModelFGY`, and `TFGY` demographic histories. They do not require a
model-specific parameter container. The two-host zoonosis study is available
separately:

```@example study
using Phylodynamics
using Phylodynamics.Zoonosis

parameters = model2_params()
result = solve_coalescent_odes(parameters; n1=20.0, n2=10.0)
(spillovers=result.summary.J_true_12, parsimony=result.summary.J_pars_1to2)
```

## Breaking namespace change

Previously, `using Phylodynamics` also imported the study's API. Now add
`using Phylodynamics.Zoonosis` for `ZoonosisParams`, `model1_params` through
`model3_params`, `build_zoonosis_*`, the study-specific trajectory/tree helpers,
`zoonosis_demographic_process`, `solve_forward_epidemic`, canonical summaries,
parameter sweeps, nested-model inference, and `plot_fig*`/`generate_all_figures`.
These names are no longer defined in the parent module.

The study extends the generic `ModelFGY`, `solve_coalescent_odes`, and
`simulate_remaster_tree` functions for `ZoonosisParams`; their reaction-system
and demographic-trajectory methods remain generic. Binary parsimony and Mk2
inference are character models, not host-specific models, and remain in the
parent module.

## Generic forward-backward simulation

`ReMASTERStepper` now executes the supplied reaction system, preserving arbitrary
species names, parameters, deme ordering, and physical migration. It no longer
reconstructs the two-host study model with fallback parameter values.
`simulate_remaster_tree(rs, infected, u0, tspan, p, samples)` also returns
all-pair `jump_counts` and `jump_times`.
Coalescent trees also preserve their model's deme order instead of inferring
reservoir/human labels, and Newick metadata writes integer states without
collapsing states beyond 2 into the first deme.

The direct SSA supports time-homogeneous propensities and integer initial
populations. Samples must be contemporaneous at `tspan[2]`; trajectories are
retried until every requested sample count is available. If multiple sampled
founders remain at `tspan[1]`, the SSA tree methods fail explicitly rather than inserting
unmodelled coalescences or assigning a reservoir root.

## Study API reference

```@docs
Phylodynamics.Zoonosis
```

```@autodocs
Modules = [Phylodynamics.Zoonosis]
Private = false
Order = [:type, :function]
```
