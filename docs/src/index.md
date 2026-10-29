# Phylodynamics.jl

Phylodynamics.jl links compartmental epidemiological models to phylogenetic
simulation and inference using ModelingToolkit.jl and Catalyst.jl reaction systems.

The package provides deterministic structured coalescent and stochastic
transmission-tree simulation, analytical coalescent ODEs, parsimony and Mk2
ancestral state reconstruction, and structured coalescent likelihoods.

## Getting started

Install from Julia's package prompt:

```julia-repl
pkg> add https://github.com/epirecipes/Phylodynamics.jl
```

The examples throughout this site run during every documentation build:

```@example quickstart
using Phylodynamics

F = [3.0 0.4 0.1; 0.2 2.0 0.3; 0.1 0.2 1.0]
G = zeros(3, 3)
Y = [40.0, 30.0, 20.0]
history = TFGY([5.0, 0.0], [F, F], [G, G], [Y, Y],
    [zeros(3), zeros(3)], ["deme_a", "deme_b", "deme_c"])
result = solve_coalescent_odes(history, [5, 4, 3];
    states=(:jumps, :fitch))
@assert length(result.summary.lineages) == 3
(lineages=result.summary.lineages, jumps=result.summary.jumps,
    parsimony=result.summary.fitch_score)
```

Follow the [tutorials](vignettes/index.md) for complete executed analyses, or consult the
[API reference](api.md) for exported types and functions.

The [multi-deme coalescent equations](coalescent-moments.md) explain selectable
analytical state groups, directed Sankoff parsimony, and symbolic verification.

The [zoonosis study](zoonosis.md) is a separate, explicitly imported submodule,
`Phylodynamics.Zoonosis`. Its parameter types, study sweeps, and publication
figures are not part of the generic namespace.

The [source repository](https://github.com/epirecipes/Phylodynamics.jl) includes
installation instructions, tests, and scripts to reproduce the study figures.
