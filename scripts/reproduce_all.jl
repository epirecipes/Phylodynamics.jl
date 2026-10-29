#!/usr/bin/env julia
# ==============================================================================
# scripts/reproduce_all.jl
# Master script for Phylodynamics.jl that runs all analytical coalescent ODE
# sweeps, ModelingToolkit reaction-network simulations, parsimony/Mk2/structured
# coalescent inference, and CairoMakie figure generation (fig1.pdf - fig9.pdf).
# ==============================================================================

using Phylodynamics
using StableRNGs
using DataFrames
using CSV
using Printf

function main()
    println("======================================================================")
    println(" Phylodynamics.jl: Complete Simulation & Inference Pipeline")
    println("======================================================================")

    rng = StableRNG(2026)
    data_dir = joinpath(@__DIR__, "..", "julia_output", "data")
    fig_dir  = joinpath(@__DIR__, "..", "julia_output", "figures")

    t0 = time()
    res = generate_all_figures(; data_dir=data_dir, fig_dir=fig_dir, rng=rng)
    elapsed = round(time() - t0, digits=2)

    println("\n--- Canonical Coalescent ODE Summary (n1 = 100, n2 = 100, T = 12.85) ---")
    show(stdout, "text/plain", res.ode_summary[:, [:model, :J_true_12, :J_true_21, :J_pars_1to2, :J_pars_2to1, :J_pars_fitch, :r_assort, :IS_norm]])
    println("\n\n--- Simulated Tree Inference Summary ---")
    show(stdout, "text/plain", res.sim_summary)
    println("\n\nAll CSV data written to: $data_dir")
    println("All 9 PDF/PNG figures written to: $fig_dir (completed in $(elapsed)s)")
end

main()
