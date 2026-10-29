# Executable tutorials

These six worked examples progress from defining reaction networks to simulation
and inference with custom multi-deme models. Every Julia cell is run as a named
Documenter `@example` block on every build. Code, tables, figures, and the complete
tutorial text appear together in the documentation; errors fail the build.

| Vignette | Topics |
|:---------|:-------|
| [1. Reaction networks](01-reaction-networks/index.md) | ModelingToolkit/Catalyst models, reaction classification, ODE and SSA trajectories |
| [2. Tree simulation](02-tree-simulation/index.md) | Structured coalescent, forward SSA, and forward-backward tree simulation |
| [3. Coalescent and parsimony ODEs](03-coalescent-parsimony-odes/index.md) | Analytical coalescent ODEs, parsimony priors, and tree topology |
| [4. Mk2 versus parsimony](04-mk2-vs-parsimony/index.md) | Ancestral state reconstruction, maximum likelihood, and stochastic mapping |
| [5. Structured coalescent inference](05-structured-coalescent-inference/index.md) | Likelihood profiles, ancestral character estimation, and island models |
| [6. Custom compartmental models](06-custom-compartmental-models/index.md) | SEIR and wildlife-livestock-human multi-deme models |

The authoritative tutorial sources are the Quarto documents in
[`vignettes/`](https://github.com/epirecipes/Phylodynamics.jl/tree/main/vignettes).
The build converts their Julia cells to `@example` blocks and executes them
against the checked-out package; it does not publish pre-rendered HTML outputs.
Each tutorial shares state between its cells but has its own execution namespace.
