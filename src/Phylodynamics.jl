"""
    Phylodynamics

Julia package for epidemiological phylodynamics, providing:
1. Structured coalescent simulation and demographic ODE integration (`ModelFGY`, `SampleConfiguration`, `SimTree`, `solveodes`, `tonewick`, `distancematrix`).
2. Structured coalescent likelihoods and ancestral character estimation (`DatedTree`, `TFGY`, `build_demographic_process`, `colik` / `colik_pik_fgy` with `PL2`, `PL1`, and `QL` approximations, `solvePikL1`, `solvePikL0`, `solveQAL`, and forward-backward `ace`).
3. Equilibrium island-model structured coalescent inference (`gen_island_model`, `tips2states`, `island_coalescent_ml`, `confint_island_coalescent`, `boot_island_coalescent`).
4. `ModelingToolkit.jl` / `Catalyst.jl` reaction-network tree generation (`TreeProblem`, `CoalescentStepper`, `ForwardSSAStepper`, `ForwardBackwardSSAStepper`), forward ODE (`OrdinaryDiffEq.jl`), and exact stochastic SSA (`JumpProcesses.jl`) simulation.
5. Analytical backward coalescent ODEs for parsimony (`1 -> 2` Camin–Sokal, `2 -> 1` reverse Camin–Sokal, and 6-state bidirectional Fitch), Mk2 continuous-time Markov chain ML & conditional stochastic character mapping, and figure generation (`CairoMakie.jl`).
"""
module Phylodynamics

using CSV
using CairoMakie
using Catalyst
using DataFrames
using Dates
using DelimitedFiles
using DiffEqBase
using DiffEqCallbacks
using Distributions
using ForwardDiff
using Interpolations
using JumpProcesses
using LinearAlgebra
using MacroTools
using ModelingToolkit
using Optim
using OrdinaryDiffEqCore
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqTsit5
using OrdinaryDiffEqVerner
using Printf
using Random
using SciMLBase
using StableRNGs
using Statistics
using StatsBase
using Symbolics

# Export core types and tree utilities
export ZoonosisParams,
       model1_params,
       model2_params,
       model3_params,
       PhyloTree,
       ntips,
       ntotal_nodes,
       root_node,
       read_newick,
       read_newick_file,
       write_newick,
       write_newick_file,
       collapse_singletons,
       postorder_edges,
       node_heights,
       node_depths,
       cherry_assortativity,
       sackin_index,

       # Coalescent.jl API
       RXN_BIRTH,
       RXN_MIG,
       RXN_DEATH,
       RXN_DYNVAR,
       SAMPLE,
       COALESCENT,
       MIGRATION,
       RECOMBINATION,
       ALGO_MARKOV,
       ALGO_STATIONARY,
       CoalescentReaction,
       ModelFGY,
       solveodes,
       SampleConfiguration,
       Event,
       SimTree,
       tonewick,
       distancematrix,
       simtree_to_phylotree,
       phylotree_to_simtree,

       # Structured Coalescent & Island Model Likelihoods
       DatedTree,
       TFGY,
       build_demographic_process,
       zoonosis_demographic_process,
       eventTimes2extant,
       update_alpha0,
       rco_finite_size_correction2!,
       solvePikL0,
       solvePikL1,
       solveQAL,
       solveQfwd0,
       colik,
       colik_pik,
       colik_pik_fgy,
       ace,
       gen_island_model,
       tips2states,
       PhylandFit,
       PhylandProfile,
       phylandml,
       confint_phylandml,
       boot_phylandml,
       island_coalescent_ml,
       confint_island_coalescent,
       boot_island_coalescent,

       # ModelingToolkit.jl / Catalyst.jl Tree Simulation
       PhyloReactionSpec,
       TreeProblem,
       CoalescentStepper,
       ForwardSSAStepper,
       ReMASTERStepper,
       ForwardBackwardSSAStepper,
       build_zoonosis_reaction_system,
       build_zoonosis_catalyst_network,
       classify_phylodynamic_reactions,
       reaction_system_to_fgy,
       simulate_ode_trajectory,
       simulate_stochastic_trajectory,
       simulate_tracked_ssa,
       simulate_forward_master_tree,
       simulate_remaster_tree,
       simulate_forward_lineage_tree,
       simulate_forward_backward_tree,
       simulate_coalescent_tree,

       # Analytical Coalescent ODEs
       solve_forward_epidemic,
       solve_coalescent_odes,
       compute_canonical_summary,
       sweep_sampling_time,
       sweep_sampling_fraction,
       sweep_R22,
       sweep_psi,

       # Inference (Parsimony, Mk2 ML + Stochastic Mapping, Structured Coalescent ML & LRT)
       parsimony_counts,
       fit_mk2_ml,
       structured_coalescent_loglik,
       fit_structured_coalescent,
       profile_structured_coalescent,

       # Figure Generation (CairoMakie.jl)
       generate_all_figures,
       plot_fig1_epidemic_and_coalescent,
       plot_fig2_simulated_trees,
       plot_fig3_sampling_time_bias,
       plot_fig4_sampling_fraction_bias,
       plot_fig5_R22_sweep,
       plot_fig6_parsimony_vs_mk2,
       plot_fig7_structured_coalescent_profiles,
       plot_fig8_topology_diagnostics,
       plot_fig9_reverse_zoonosis_sweep

include("types_and_trees.jl")
include("coalescent_base.jl")
include("phydynr_phyland.jl")
include("reaction_network.jl")
include("coalescent_odes.jl")
include("inference.jl")
include("figures.jl")

end # module Phylodynamics
