"""
    Phylodynamics.Zoonosis

Two-host zoonosis study models, legacy 37-state equations, parameter sweeps,
nested-model inference, and publication figures. Import explicitly with
`using Phylodynamics.Zoonosis`; generic reaction-system and demographic-trajectory
APIs remain in `Phylodynamics` and do not require `ZoonosisParams`.
"""
module Zoonosis

using ..Phylodynamics
import ..Phylodynamics: ModelFGY, solve_coalescent_odes, simulate_remaster_tree
using ..Phylodynamics: _tree_layout_coords
using CSV
using CairoMakie
using Catalyst
using DataFrames
using Distributions
using ModelingToolkit
using Optim
using OrdinaryDiffEqTsit5
using Printf
using Random
using SciMLBase
using StableRNGs
using Statistics
using StatsBase
using JumpProcesses

export ZoonosisParams, model1_params, model2_params, model3_params,
       build_zoonosis_reaction_system, build_zoonosis_catalyst_network,
       simulate_ode_trajectory, simulate_stochastic_trajectory,
       simulate_coalescent_tree, solve_forward_epidemic,
       zoonosis_demographic_process, compute_canonical_summary,
       sweep_sampling_time, sweep_sampling_fraction, sweep_R22, sweep_psi,
       structured_coalescent_loglik, fit_structured_coalescent,
       profile_structured_coalescent, generate_all_figures,
       plot_fig1_epidemic_and_coalescent, plot_fig2_simulated_trees,
       plot_fig3_sampling_time_bias, plot_fig4_sampling_fraction_bias,
       plot_fig5_R22_sweep, plot_fig6_parsimony_vs_mk2,
       plot_fig7_structured_coalescent_profiles, plot_fig8_topology_diagnostics,
       plot_fig9_reverse_zoonosis_sweep

include("parameters.jl")
include("reaction_network.jl")
include("demography.jl")
include("coalescent_odes.jl")
include("inference.jl")
include("figures.jl")

end
