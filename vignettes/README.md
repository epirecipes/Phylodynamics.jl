# Vignettes

GitHub Pages publishes these sources as full Documenter tutorials. Every Julia
cell is converted to a shared, named `@example` block and executed on each build;
errors fail the build. To build the tutorials together with the API reference:

```sh
cd Phylodynamics.jl
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

Each subdirectory holds one Quarto document (Julia engine) that can be rendered on its own:

```sh
cd Phylodynamics.jl/vignettes
julia --project=. -e 'using Pkg; Pkg.develop(path=".."); Pkg.instantiate()'
quarto render 01-reaction-networks/reaction-networks.qmd   # one vignette
quarto render                                              # all of them
```

The vignettes progress from basic reaction-network specification and epidemic simulation to advanced structured coalescent inference and custom multi-deme compartmental models:

| # | Vignette | What it shows |
|:--|:---------|:--------------|
| 1 | `01-reaction-networks` | Defining two-host zoonotic models with `ModelingToolkit.jl` and `Catalyst.jl`; automatic stoichiometric classification of transmission ($F_{kl}$), migration ($G_{kl}$), removal ($D_k$), and demographic reactions; deterministic ODE (`Tsit5`) and exact Gillespie SSA (`JumpProcesses.jl`) trajectories across Models 1, 2, and 3 |
| 2 | `02-tree-simulation` | Formulating a `TreeProblem` with `Pair` maps (`[Y1 => 80, Y2 => 40]`) and simulating phylogenies via three engines: deterministic structured coalescent (`CoalescentStepper`), forward-in-time individual-lineage SSA (`ForwardSSAStepper`), and forward-backward reaction event tracing (`ForwardBackwardSSAStepper`) |
| 3 | `03-coalescent-parsimony-odes` | The 37-state backward structured coalescent ODEs (`solve_coalescent_odes`): expected true jumps ($J^*_{12}, J^*_{21}$), the three parsimony flow priors ($1\to 2$ Camin–Sokal, $2\to 1$ reverse Camin–Sokal, and 6-state bidirectional Fitch), cherry assortativity $r$, normalized Sackin tree imbalance $I_S^*$, and topological detection of reverse zoonosis ($\Delta J = J^{(0)} - J_F > 0$) |
| 4 | `04-mk2-vs-parsimony` | Ancestral state reconstruction on phylogenies: comparing the three parsimony priors (`parsimony_counts`) against 1-rate unidirectional and 2-rate bidirectional continuous-time Markov ("mugration" / Mk2) maximum likelihood (`fit_mk2_ml`), conditional stochastic character mapping, and opposite sample-size scaling of $\hat{q}_{12}$ versus $\hat{q}_{21}$ |
| 5 | `05-structured-coalescent-inference` | Structured coalescent inference and ancestral character estimation: building a demographic process directly from a `ReactionSystem` (`build_demographic_process`), `DatedTree`, structured coalescent likelihoods (`colik` with `:PL2`, `:PL1`, and `:QL`), 1D/2D profile likelihoods, LRTs for human-to-human spread and reverse zoonosis, forward-backward ancestral character estimation (`ace`), and equilibrium island-model inference (`island_coalescent_ml`, `confint_island_coalescent`, `boot_island_coalescent`) |
| 6 | `06-custom-compartmental-models` | Generalizing beyond two-host SIRS: building custom multi-stage SEIR (`[E, I]`) and 3-deme wildlife–livestock–human bridge-reservoir (`[Y1, Y2, Y3]`) `ReactionSystem` models, automatic stoichiometric `TFGY` construction, multi-deme genealogy simulation, and 3-state ancestral host reconstruction |
