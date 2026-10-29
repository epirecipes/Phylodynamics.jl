"""
    Phylodynamics

Julia package for epidemiological phylodynamics, providing:
1. Structured coalescent simulation and demographic ODE integration (`ModelFGY`, `SampleConfiguration`, `SimTree`, `solveodes`, `tonewick`, `distancematrix`).
2. Structured coalescent likelihoods and ancestral character estimation (`DatedTree`, `TFGY`, `build_demographic_process`, `colik` / `colik_pik_fgy` with `PL2`, `PL1`, and `QL` approximations, `solvePikL1`, `solvePikL0`, `solveQAL`, and forward-backward `ace`).
3. Equilibrium island-model structured coalescent inference (`gen_island_model`, `tips2states`, `island_coalescent_ml`, `confint_island_coalescent`, `boot_island_coalescent`).
4. `ModelingToolkit.jl` / `Catalyst.jl` reaction-network tree generation (`TreeProblem`, `CoalescentStepper`, `ForwardSSAStepper`, `ForwardBackwardSSAStepper`), forward ODE (`OrdinaryDiffEq.jl`), and exact stochastic SSA (`JumpProcesses.jl`) simulation.
5. Selectable multi-deme analytical coalescent moments, Fitch and directed/weighted Sankoff parsimony, and binary Mk2 character inference. The two-host study and its reproduction utilities are available through an explicit `using Phylodynamics.Zoonosis`.
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
export PhyloTree,
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
       classify_phylodynamic_reactions,
       reaction_system_to_fgy,
       simulate_tracked_ssa,
       simulate_forward_master_tree,
       simulate_remaster_tree,
       simulate_forward_lineage_tree,
       simulate_forward_backward_tree,

       # Analytical Coalescent ODEs
       solve_coalescent_odes,

       # Binary character inference
       parsimony_counts,
       fit_mk2_ml

include("types_and_trees.jl")
include("coalescent_base.jl")
include("phydynr_phyland.jl")
include("reaction_network.jl")
include("remaster.jl")
include("coalescent_moments.jl")
include("inference.jl")
include("tree_layout.jl")
include("zoonosis/Zoonosis.jl")

end # module Phylodynamics
