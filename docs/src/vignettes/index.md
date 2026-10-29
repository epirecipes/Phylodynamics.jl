# Vignettes

These six worked examples progress from defining reaction networks to simulation
and inference with custom multi-deme models. Each link opens a standalone rendered
Quarto vignette, including its code, results, and figures.

| Vignette | Topics |
|:---------|:-------|
| [1. Reaction networks](01-reaction-networks/reaction-networks.html) | ModelingToolkit/Catalyst models, reaction classification, ODE and SSA trajectories |
| [2. Tree simulation](02-tree-simulation/tree-simulation.html) | Structured coalescent, forward SSA, and forward-backward tree simulation |
| [3. Coalescent and parsimony ODEs](03-coalescent-parsimony-odes/coalescent-parsimony-odes.html) | Analytical coalescent ODEs, parsimony priors, and tree topology |
| [4. Mk2 versus parsimony](04-mk2-vs-parsimony/mk2-vs-parsimony.html) | Ancestral state reconstruction, maximum likelihood, and stochastic mapping |
| [5. Structured coalescent inference](05-structured-coalescent-inference/structured-coalescent-inference.html) | Likelihood profiles, ancestral character estimation, and island models |
| [6. Custom compartmental models](06-custom-compartmental-models/custom-compartmental-models.html) | SEIR and wildlife-livestock-human multi-deme models |

The rendered vignettes are published from the HTML files committed in
[`vignettes/`](https://github.com/epirecipes/Phylodynamics.jl/tree/main/vignettes).
After changing a vignette, render it with Quarto and commit both the `.qmd` source
and the updated `.html` before publishing.
