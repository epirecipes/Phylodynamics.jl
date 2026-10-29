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
using StableRNGs

tree = simulate_coalescent_tree(model2_params(); n1=20, n2=10, rng=StableRNG(42))
counts = parsimony_counts(tree)
@assert ntips(tree) == 30
(tips=ntips(tree), directional_changes=counts.J_pars_1to2,
    fitch_changes=counts.J_pars_fitch)
```

Follow the [tutorials](vignettes/index.md) for complete executed analyses, or consult the
[API reference](api.md) for exported types and functions.

The [multi-deme coalescent equations](coalescent-moments.md) explain selectable
analytical state groups, directed Sankoff parsimony, and symbolic verification.

The [source repository](https://github.com/epirecipes/Phylodynamics.jl) includes
installation instructions, tests, and scripts to reproduce the study figures.
