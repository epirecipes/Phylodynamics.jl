# Phylodynamics.jl

Phylodynamics.jl links compartmental epidemiological models to phylogenetic
simulation and inference using ModelingToolkit.jl and Catalyst.jl reaction systems.

The package provides deterministic structured coalescent and stochastic
transmission-tree simulation, analytical coalescent ODEs, parsimony and Mk2
ancestral state reconstruction, and structured coalescent likelihoods.

## Getting started

```julia
using Pkg
Pkg.add(url="https://github.com/epirecipes/Phylodynamics.jl")
```

Follow the [vignettes](vignettes/index.md) for worked examples, or consult the
[API reference](api.md) for exported types and functions.

The [source repository](https://github.com/epirecipes/Phylodynamics.jl) includes
installation instructions, tests, and scripts to reproduce the study figures.
