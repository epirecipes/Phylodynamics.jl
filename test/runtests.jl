using Test
using Phylodynamics
using Phylodynamics.Zoonosis
using ModelingToolkit
using Catalyst
using SciMLBase
using StableRNGs
using DataFrames
using CSV

include("documentation.jl")
include("generic_framework.jl")
include("coalescent_moments.jl")
include("coalescent_symbolics.jl")

@testset "Phylodynamics.jl Complete Test Suite" begin

    @testset "1. ModelingToolkit ReactionSystem & Stoichiometric Classification" begin
        bundle = build_zoonosis_reaction_system()
        specs = classify_phylodynamic_reactions(bundle.rs, bundle.infected)
        @test length(specs) == 12
        roles = [s.role for s in specs]
        @test count(==(:transmission), roles) == 4
        @test count(==(:removal), roles) == 4
        @test count(==(:demographic), roles) == 4

        # Check parent -> child demes for the 4 transmission reactions
        trans_pairs = [(s.parent_deme, s.child_deme) for s in specs if s.role === :transmission]
        @test (1, 1) in trans_pairs # Reservoir -> Reservoir
        @test (1, 2) in trans_pairs # Spillover 1 -> 2
        @test (2, 2) in trans_pairs # Human -> Human 2 -> 2
        @test (2, 1) in trans_pairs # Reverse zoonosis 2 -> 1

        # Verify Catalyst @reaction_network compatibility as well
        rn_cat = build_zoonosis_catalyst_network()
        sp_cat = Catalyst.species(rn_cat)
        y1_cat = sp_cat[findfirst(s -> occursin("Y1", string(s)), sp_cat)]
        y2_cat = sp_cat[findfirst(s -> occursin("Y2", string(s)), sp_cat)]
        specs_cat = classify_phylodynamic_reactions(rn_cat, [y1_cat, y2_cat])
        @test count(s -> s.role === :transmission, specs_cat) == 4
    end

    @testset "2. Forward ODE & Stochastic JumpProcesses Simulation from MTK" begin
        zp1 = model1_params()
        ode_res = simulate_ode_trajectory(zp1; saveat=0.0:1.0:12.85)
        @test ode_res.sol.retcode == ReturnCode.Success
        rng = StableRNG(123)
        jmp_res = simulate_stochastic_trajectory(zp1; rng=rng)
        @test jmp_res.sol.t[end] ≈ 12.85
    end

    @testset "3. Tree Simulation: Coalescent, Forward SSA, and Forward-Backward SSA from MTK" begin
        bundle = build_zoonosis_reaction_system()
        sp = bundle.species
        pr = bundle.params
        zp2 = model2_params()
        u0 = [sp.X1 => zp2.X1_0, sp.Y1 => zp2.Y1_0, sp.X2 => zp2.X2_0, sp.Y2 => zp2.Y2_0]
        p = [
            pr.nu1 => zp2.nu1, pr.mu1 => zp2.mu1, pr.gamma1 => zp2.gamma1, pr.m1 => zp2.m1,
            pr.beta => zp2.beta, pr.nu2 => zp2.nu2, pr.mu2 => zp2.mu2, pr.gamma2 => zp2.gamma2,
            pr.m2 => zp2.m2, pr.epsilon => zp2.epsilon, pr.phi => zp2.phi, pr.psi => zp2.psi
        ]
        tspan = (0.0, 12.85)
        samples = [sp.Y1 => 30, sp.Y2 => 30]
        tprob = TreeProblem(bundle.rs, bundle.infected, u0, tspan, p, samples)

        # (a) Deterministic CoalescentStepper
        tr_coal = solve(tprob, CoalescentStepper(); rng=StableRNG(101))
        @test ntips(tr_coal) == 60
        @test tr_coal.nnode == 59
        @test tr_coal.true_jumps_12 >= 1

        # (b) Forward-in-time stochastic simulation
        tr_fwd = solve(tprob, ForwardSSAStepper(); rng=StableRNG(102))
        @test ntips(tr_fwd) == 60
        @test tr_fwd.nnode == 59
        @test tr_fwd.true_jumps_12 >= 1

        # (c) Forward-backward stochastic reaction tree simulation
        tr_rem = solve(tprob, ForwardBackwardSSAStepper(); rng=StableRNG(103))
        @test ntips(tr_rem) == 60
        @test tr_rem.nnode == 59
        @test tr_rem.true_jumps_12 >= 1

        # Round-trip Newick serialization
        nwk = write_newick(tr_rem; include_metadata=true)
        tr_parsed = read_newick(nwk)
        @test ntips(tr_parsed) == 60
    end

    @testset "4. Analytical Backward Coalescent ODEs (Models 1, 2, 3)" begin
        # Cross-check against R reference (n1 = 100, n2 = 10, T = 12.85)
        df_10 = compute_canonical_summary(; n1=100.0, n2=10.0, T_sample=12.85)
        @test nrow(df_10) == 3
        @test isapprox(df_10.J_true_12[1], 10.0; atol=0.05)
        @test isapprox(df_10.J_pars_1to2[1], 9.69; atol=0.05)
        @test isapprox(df_10.J_pars_2to1[1], 32.65; atol=0.15)
        @test isapprox(df_10.J_pars_fitch[1], 9.67; atol=0.05)

        @test isapprox(df_10.J_true_12[2], 7.12; atol=0.05)
        @test isapprox(df_10.J_pars_1to2[2], 6.96; atol=0.05)
        @test isapprox(df_10.r_assort[2], 0.488; atol=0.01)

        @test isapprox(df_10.J_true_12[3], 8.09; atol=0.05)
        @test isapprox(df_10.J_true_21[3], 4.06; atol=0.05)
        @test df_10.J_pars_2to1[3] > df_10.J_pars_1to2[3] > df_10.J_pars_fitch[3]

        # Balanced sampling (n1 = 100, n2 = 100, T = 12.85)
        df_sum = compute_canonical_summary(; n1=100.0, n2=100.0, T_sample=12.85)
        @test isapprox(df_sum.J_true_12[1], 100.0; atol=0.1)
        @test isapprox(df_sum.J_true_21[1], 0.0; atol=1e-6)
        @test df_sum.r_assort[1] < 0.0
        @test isapprox(df_sum.J_true_12[2], 28.86; atol=0.2)
        @test df_sum.r_assort[2] > 0.35
        @test df_sum.J_true_21[3] > 4.0
        @test df_sum.J_pars_2to1[3] > df_sum.J_pars_1to2[3] > df_sum.J_pars_fitch[3]
    end

    @testset "5. Parsimony, Mk2 ML, Structured Coalescent (colik/ace), and Island Coalescent" begin
        rng = StableRNG(2026)
        t1 = simulate_forward_backward_tree(model1_params(); n1=40, n2=40, rng=rng).binary_tree
        t2 = simulate_forward_backward_tree(model2_params(); n1=40, n2=40, rng=rng).binary_tree
        t3 = simulate_forward_backward_tree(model3_params(); n1=40, n2=40, rng=rng).binary_tree

        pc1 = parsimony_counts(t1)
        pc3 = parsimony_counts(t3)
        @test pc1.J_pars_1to2 >= pc1.J_pars_fitch
        @test pc3.J_fitch_12 + pc3.J_fitch_21 == pc3.J_pars_fitch

        mk3 = fit_mk2_ml(t3; n_maps=20, rng=rng)
        @test mk3.J_ml_bidir_tot >= pc3.J_pars_fitch * 0.8
        @test isfinite(mk3.loglik_unidir) && isfinite(mk3.loglik_bidir)

        # Structured coalescent colik (PL2, PL1, QL) and ace
        bdt2 = DatedTree(t2; demes=["Y1", "Y2"])
        dm2 = zoonosis_demographic_process(model2_params())
        theta2 = Dict("epsilon" => 0.01, "phi" => 4.0, "psi" => 0.0)
        ll_ql  = colik(bdt2, theta2, dm2; res=200, likelihood=:QL, timeOfOriginBoundaryCondition=false)
        ll_pl1 = colik(bdt2, theta2, dm2; res=200, likelihood=:PL1, timeOfOriginBoundaryCondition=false)
        ll_pl2 = colik(bdt2, theta2, dm2; res=200, likelihood=:PL2, timeOfOriginBoundaryCondition=false)
        @test isfinite(ll_ql) && isfinite(ll_pl1) && isfinite(ll_pl2)

        bdt2_ace = ace(bdt2, theta2, dm2; res=100, likelihood=:PL2, timeOfOriginBoundaryCondition=false)
        @test size(bdt2_ace.acestates) == (2, bdt2.n + bdt2.Nnode)
        @test all(!isnan, bdt2_ace.acestates[:, bdt2.n + 1])

        # Equilibrium island model ML
        fit_pl = island_coalescent_ml(t2; maxiters=80, likelihood=:QL)
        @test isfinite(fit_pl.loglik)
        @test nrow(fit_pl.ace) == bdt2.n + bdt2.Nnode
    end

    @testset "6. Reference Data Verification" begin
        ref_dir = isdir(joinpath(@__DIR__, "..", "zoonosis_phylodynamics", "data")) ?
            joinpath(@__DIR__, "..", "zoonosis_phylodynamics", "data") :
            joinpath(@__DIR__, "..", "..", "zoonosis_phylodynamics", "data")
        if isdir(ref_dir)
            # (a) Forward ODE trajectories (7-digit agreement)
            df_fwd_r = CSV.read(joinpath(ref_dir, "ode_forward_trajectories.csv"), DataFrame)
            for (m_r, zp) in [("Model1_Spillover", model1_params()), ("Model2_H2H", model2_params()), ("Model3_RevZoonosis", model3_params())]
                sub_r = filter(r -> r.model == m_r && r.time <= 12.85, df_fwd_r)
                sol_jl = solve_forward_epidemic(zp; dt=0.01).sol
                max_rel_Y1 = maximum(abs(sol_jl(r.time)[2] - r.Y1) / max(r.Y1, 1.0) for r in eachrow(sub_r))
                max_rel_Y2 = maximum(abs(sol_jl(r.time)[4] - r.Y2) / max(r.Y2, 1.0) for r in eachrow(sub_r))
                @test max_rel_Y1 < 1e-6
                @test max_rel_Y2 < 1e-6
            end

            # (b) Unchanged moments retain their R references. Historical directional
            # columns violate class conservation; corrected equations and trajectories
            # are checked in coalescent_symbolics.jl and coalescent_moments.jl.
            df_curve_r = CSV.read(joinpath(ref_dir, "coalescent_model2_n1_curve.csv"), DataFrame)
            for n1_val in [10, 50, 100]
                row_r = filter(r -> r.n1 == n1_val, df_curve_r)[1, :]
                s2 = solve_coalescent_odes(model2_params(); n1=Float64(n1_val), n2=10.0).summary
                s3 = solve_coalescent_odes(model3_params(); n1=Float64(n1_val), n2=10.0).summary
                @test isapprox(s2.r_assort, row_r.m2_r; atol=1e-3)
                @test isapprox(s3.J_pars_fitch, row_r.m3_J_Fitch; atol=5e-3)
                @test isapprox(s3.J_true_21, row_r.m3_Js21; atol=5e-3)
            end

            # (c) Exact parsimony (81/81 trees) and Mk2 ML agreement on ReMASTER NEXUS trees
            df_mk_r = filter(r -> r.source == "ReMASTER", CSV.read(joinpath(ref_dir, "parsimony_mk2_summary.csv"), DataFrame))
            for (m_r, prefix) in [
                ("Model1_Spillover", "remaster_Model1_Spillover"),
                ("Model2_H2H", "remaster_Model2_H2H"),
                ("Model3_RevZoonosis", "remaster_Model3_RevZoonosis")
            ]
                trees = read_newick_file(joinpath(ref_dir, "$(prefix)_n1_100_n2_10.binary.trees"))
                @test length(trees) == 9
                for (tid, tr_raw) in enumerate(trees)
                    tr = collapse_singletons(tr_raw)
                    pc = parsimony_counts(tr)
                    mk = fit_mk2_ml(tr; n_maps=5)
                    row_r = filter(r -> r.model == m_r && r.n1 == 100 && r.tree_id == tid, df_mk_r)[1, :]
                    @test pc.J_pars_1to2 == row_r.cs_jumps
                    @test pc.J_pars_2to1 == row_r.cs21_jumps
                    @test pc.J_pars_fitch == row_r.fitch_jumps
                    @test isapprox(sum(tr.edge_length), row_r.tree_length; atol=1e-9)
                    @test isapprox(log(mk.q_unidir), row_r.log_q12_mle; atol=1e-4)
                    @test isapprox(mk.loglik_unidir, row_r.mk2_loglik; atol=1e-6)
                end
            end
        end
    end

end
