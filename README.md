# Phylodynamics.jl

A Julia package for forward and backward phylodynamic simulation, analytical structured coalescent ODEs, phylogeographic ancestral state reconstruction, and maximum-likelihood phylodynamic inference, built directly on **ModelingToolkit.jl** and **Catalyst.jl** reaction systems.

**Documentation:** [Vignettes and API reference](https://epirecip.es/Phylodynamics.jl/).

The core API is model-independent: supply a reaction system, its infected
species, initial conditions, parameter map, and samples, or a `TFGY` demographic
history. No `ZoonosisParams` is required. The original two-host study is now
separate: add **`using Phylodynamics.Zoonosis`** for its presets, sweeps,
study-specific inference, and publication figures. Those names are no longer
available from `using Phylodynamics` alone; see the
[migration guide](https://epirecip.es/Phylodynamics.jl/zoonosis/).

---

## Key Features

1. **SciML / ModelingToolkit.jl Reaction System Syntax**
   - Define epidemiological models using standard `ModelingToolkit.jl` (`@independent_variables`, `@parameters`, `@species`, `Reaction`, `ReactionSystem`) or `Catalyst.jl` (`@reaction_network`) syntax, together with `[species => value]` `Pair` maps for initial conditions `u0`, parameters `p`, and tip sampling configurations `samples`.
   - `classify_phylodynamic_reactions(rs, infected_species)` automatically inspects reaction stoichiometry on the infected compartments (`[Y1, Y2, ...]`) to classify reactions into:
     - `:transmission` ($F_{kl}$: parent deme $k \to$ child deme $l$, net infected $+1$),
     - `:migration` ($G_{kl}$: deme $k \to$ deme $l$, net infected $0$),
     - `:removal` ($D_k$: deme $k \to \emptyset$, net infected $-1$), and
     - `:demographic` (uninfected compartment turnover).

2. **Three Unified Tree Simulation Engines via `TreeProblem(rs, infected, u0, tspan, p, samples)`**
   - **`CoalescentStepper()`** (*Deterministic backward-in-time structured coalescent*): Converts `ReactionSystem` into an `ODESystem`, solves the forward demographic trajectory with `OrdinaryDiffEqTsit5`, and simulates genealogies backward in time via cumulative-hazard `ContinuousCallback` integration (`ModelFGY`, `SampleConfiguration`, `SimTree`).
   - **`ForwardSSAStepper()`** (*Forward-in-time individual-lineage stochastic simulation*): Runs an exact Gillespie stochastic simulation while tracking individual infected host lineages forward in time (`X1 + Y1 --> 2Y1`, `X2 + Y1 --> Y2 + Y1`, etc.), building the true transmission genealogy forward in time and pruning to sampled tips.
   - **`ForwardBackwardSSAStepper()`** (*Forward-backward stochastic reaction tree simulation*): Runs an exact Gillespie SSA from `Catalyst.jumpratelaw`s and `netstoichmat(rs)`, then reconstructs both the full singleton-annotated reaction tree (`[&type="Y1",time=...]`) and the collapsed binary phylogeny backward in time.

3. **Structured Coalescent Likelihoods, Ancestral Character Estimation, and Island Models**
   - **Demographic & Genealogical Models**: `ModelFGY`, `SampleConfiguration`, `CoalescentReaction`, `solveodes`, `Event`, `SimTree`, `tonewick`, `distancematrix`.
   - **Time-Varying Structured Coalescent Inference**: `DatedTree`, `TFGY`, `build_demographic_process`, `colik` / `colik_pik` / `colik_pik_fgy` (`:PL2`, `:PL1`, and `:QL` approximations via `solvePikL1`, `solvePikL0`, `solveQAL`), and forward-backward ancestral character estimation `ace` (`solveQfwd0`).
   - **Equilibrium Island Coalescent Inference**: `gen_island_model`, `tips2states`, `island_coalescent_ml`, `confint_island_coalescent`, and `boot_island_coalescent`.

4. **Selectable Multi-Deme Coalescent & Parsimony ODEs (`solve_coalescent_odes`)**
   - Accepts a `TFGY` trajectory or a Catalyst/ModelingToolkit reaction system, with optional lineage, directional-jump, coalescence, branch-length, cherry, Sackin, Fitch, and directed/weighted Sankoff state groups. Parsimony closure is opt-in and bounded by explicit state budgets; observed character labels can group multiple demographic demes.
   - `Phylodynamics.Zoonosis` supplies the two-host, 37-state specialization for Models 1–3. Symbolics.jl regenerates all 37 derivatives from the general event equations, and a local-venv SymPy verifier independently checks the algebra and small-tree Sankoff scores.
   - **Correction:** the legacy directional Camin-Sokal class equations previously failed to conserve lineage counts. Those terms are fixed; old directional-parsimony summaries need regeneration. See the [derivation and verification guide](https://epirecip.es/Phylodynamics.jl/coalescent-moments/).

5. **Publication Figures via `CairoMakie.jl`**
   - After `using Phylodynamics.Zoonosis`, `generate_all_figures()` reproduces the study's summary figures (`fig1`–`fig9` in `.pdf` and `.png`) and CSV tables in `julia_output/`.

---

## Quick Start: ModelingToolkit Reaction System to Phylogeny & Inference

```julia
using Phylodynamics
using ModelingToolkit
using Catalyst
using OrdinaryDiffEqTsit5
using JumpProcesses
using StableRNGs

# 1. Define the two-host zoonosis reaction network with ModelingToolkit / Catalyst
t = ModelingToolkit.t_nounits
@parameters nu1 mu1 gamma1 m1 beta nu2 mu2 gamma2 m2 epsilon phi psi
@species X1(t) Y1(t) X2(t) Y2(t)

rxs = [
    Reaction(nu1, nothing, [X1], nothing, [1]),
    Reaction(mu1, [X1], nothing, [1], nothing),
    Reaction(beta, [X1, Y1], [Y1], [1, 1], [2]),
    Reaction(gamma1, [Y1], nothing, [1], nothing),
    Reaction(mu1 + m1, [Y1], nothing, [1], nothing),
    Reaction(nu2, nothing, [X2], nothing, [1]),
    Reaction(mu2, [X2], nothing, [1], nothing),
    Reaction(epsilon * beta, [X2, Y1], [Y2, Y1], [1, 1], [1, 1]),
    Reaction(phi * beta, [X2, Y2], [Y2], [1, 1], [2]),
    Reaction(psi * beta, [X1, Y2], [Y1, Y2], [1, 1], [1, 1]),
    Reaction(gamma2, [Y2], nothing, [1], nothing),
    Reaction(mu2 + m2, [Y2], nothing, [1], nothing)
]
@named rs = ReactionSystem(rxs, t, [X1, Y1, X2, Y2],
    [nu1, mu1, gamma1, m1, beta, nu2, mu2, gamma2, m2, epsilon, phi, psi];
    combinatoric_ratelaws = false)
rs = complete(rs)

# 2. Specify initial conditions, parameters (Model 3), and tip samples using Pair maps
u0 = [X1 => 9999.0, Y1 => 1.0, X2 => 1000.0, Y2 => 0.0]
p  = [nu1 => 2.74, mu1 => 2.74e-4, gamma1 => 0.1, m1 => 0.0,
      beta => 1.00274e-4, nu2 => 1000/(70*365), mu2 => 1/(70*365),
      gamma2 => 0.01, m2 => 0.1, epsilon => 0.01, phi => 4.0, psi => 4.0]
tspan   = (0.0, 12.85)
samples = [Y1 => 100, Y2 => 30]

# 3. Formulate a TreeProblem and solve with any of the three engines:
tprob = TreeProblem(rs, [Y1, Y2], u0, tspan, p, samples)

# (a) Deterministic structured coalescent
tree_coal = solve(tprob, CoalescentStepper(); rng = StableRNG(42))

# (b) Forward-in-time individual-lineage Gillespie SSA
tree_fwd  = solve(tprob, ForwardSSAStepper(); rng = StableRNG(42))

# (c) Forward-backward stochastic reaction tree
tree_fb   = solve(tprob, ForwardBackwardSSAStepper(); rng = StableRNG(42))

# 4. Run parsimony, Mk2 ML + stochastic mapping, and structured coalescent ML
pars = parsimony_counts(tree_coal)
mk2  = fit_mk2_ml(tree_coal; rng = StableRNG(42))
dm   = build_demographic_process(rs, [Y1, Y2], u0, tspan, p)
bdt  = DatedTree(tree_coal; demes = ["Y1", "Y2"])
ll   = colik(bdt, Dict("epsilon" => 0.01, "phi" => 4.0, "psi" => 4.0), dm;
             t0 = 0.0, res = 400, timeOfOriginBoundaryCondition = false, likelihood = :PL2)
```

---

## Running Tests & Full Simulation Pipeline

To run the test suite:
```bash
julia --project=. test/runtests.jl
```

To generate all analytical coalescent ODE tables, simulated phylogenies, inference summaries, and figures (`julia_output/data/` and `julia_output/figures/fig1.pdf`–`fig9.pdf`):
```bash
julia --project=. scripts/reproduce_all.jl
```

---

## Quarto Vignettes (`vignettes/`)

Six progressive Quarto sources are provided in `vignettes/`. The documentation
build converts every Julia cell into a Documenter `@example` block, executes all
six tutorials, and publishes their code, results, figures, and the complete API
reference together. A failing example fails the build; committed standalone HTML
is not used for GitHub Pages.

| # | Vignette | Description |
|:--|:---------|:------------|
| 1 | `01-reaction-networks/reaction-networks.qmd` | Defining two-host zoonotic models in `ModelingToolkit.jl` and `Catalyst.jl`, automatic stoichiometric classification (`classify_phylodynamic_reactions`), deterministic ODEs (`Tsit5`), and exact Gillespie SSA (`JumpProcesses.jl`) trajectories |
| 2 | `02-tree-simulation/tree-simulation.qmd` | Simulating phylogenies via `TreeProblem` across all three engines: deterministic structured coalescent (`CoalescentStepper`), forward-in-time individual-lineage SSA (`ForwardSSAStepper`), and forward-backward reaction event tracing (`ForwardBackwardSSAStepper`) |
| 3 | `03-coalescent-parsimony-odes/coalescent-parsimony-odes.qmd` | Analytical 37-state backward structured coalescent ODEs (`solve_coalescent_odes`), the three parsimony flow priors ($1\to 2$ CS, $2\to 1$ CS, and 6-state Fitch), cherry assortativity $r$, normalized Sackin index $I_S^*$, and reverse-zoonosis topological excess $\Delta J = J^{(0)} - J_F > 0$ |
| 4 | `04-mk2-vs-parsimony/mk2-vs-parsimony.qmd` | Comparing parsimony (`parsimony_counts`) against 1-rate unidirectional and 2-rate bidirectional Mk2 maximum likelihood (`fit_mk2_ml`), opposite sample-size scaling of $\hat{q}_{12}$ vs. $\hat{q}_{21}$, and stochastic mapping jump-time distributions |
| 5 | `05-structured-coalescent-inference/structured-coalescent-inference.qmd` | Structured coalescent likelihoods and ancestral character estimation: `build_demographic_process`, `DatedTree`, `colik` (`:PL2`, `:PL1`, `:QL`), 1D/2D profile likelihoods, LRTs for $\phi$ and $\psi$, forward-backward `ace` internal-node posteriors, and `island_coalescent_ml` equilibrium island-model inference |
| 6 | `06-custom-compartmental-models/custom-compartmental-models.qmd` | Custom multi-deme compartmental models beyond two-host SIRS: a 2-deme Exposed–Infectious SEIR model (`[E, I]`) and a 3-deme wildlife–livestock–human bridge-reservoir chain (`[Y1, Y2, Y3]`) with 3-state `ace` reconstruction |

To build the complete executable documentation from the repository root:
```bash
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

Alternatively, render the original Quarto sources:
```bash
cd vignettes
quarto render 01-reaction-networks/reaction-networks.qmd
quarto render
```
