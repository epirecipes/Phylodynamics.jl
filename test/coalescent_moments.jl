using Test
using Phylodynamics
using Phylodynamics.Zoonosis
using LinearAlgebra
using Catalyst
using DataFrames

function constant_moment_tfgy(F, G, Y; duration=2.0)
    m = length(Y)
    return TFGY([duration, 0.0], [Float64.(F), Float64.(F)],
        [Float64.(G), Float64.(G)], [Float64.(Y), Float64.(Y)],
        [zeros(m), zeros(m)], ["Y$i" for i in 1:m])
end

@testset "Selectable multi-deme coalescent moments" begin
    @testset "One-deme analytic solution and root stopping" begin
        tfgy = constant_moment_tfgy(fill(10.0, 1, 1), zeros(1, 1), [10.0]; duration=10.0)
        minimal = solve_coalescent_odes(tfgy, [4.0]; dt=0.25)
        @test minimal.layout.nstates == 1
        @test Set(keys(minimal.layout.groups)) == Set([:lineages])
        @test minimal.summary.reached_mrca
        @test last(minimal.sol.t) ≈ 7.5 atol=1e-7
        @test minimal.summary.lineages ≈ [1.0] atol=1e-9
        @test all(isapprox(u[1], 4 / (1 + 0.4s); atol=1e-8)
            for (s, u) in zip(minimal.sol.t, minimal.sol.u))
        full = solve_coalescent_odes(tfgy, [4.0];
            states=(:coalescences, :branch_lengths, :cherries, :sackin, :fitch), dt=0.25)
        @test full.summary.coalescences[1, 1] ≈ 3.0 atol=1e-7
        @test only(full.summary.branch_lengths) ≈ 10log(4) atol=1e-7
        @test only(full.summary.descendants) ≈ 4.0
        @test full.summary.sackin ≈ 8log(4) atol=1e-7
        @test only(full.summary.cherries) ≈ 4 / 3 * (1 - 4.0^-3) atol=1e-7
        @test full.summary.fitch_score == 0
    end

    @testset "Migration-only transport" begin
        F, G = zeros(2, 2), [0.0 2.0; 0.0 0.0]
        tfgy = constant_moment_tfgy(F, G, [10.0, 20.0])
        result = solve_coalescent_odes(tfgy, [0.0, 4.0];
            states=(:jumps, :coalescences, :sackin, :cherries, :fitch))
        moved = 4 * (1 - exp(-0.2))
        @test result.summary.lineages ≈ [moved, 4 - moved] atol=1e-9
        @test result.summary.descendants ≈ result.summary.lineages atol=1e-9
        @test result.summary.jumps[1, 2] ≈ moved atol=1e-9
        @test result.summary.coalescences == zeros(2, 2)
        @test result.summary.cherries == zeros(2, 2)
        @test result.summary.sackin == result.summary.fitch_score == 0
        @test !result.summary.reached_mrca
    end

    @testset "Nonuniform demographic interpolation and restricted interval" begin
        tfgy = TFGY([5.0, 4.0, 0.0], [fill(f, 1, 1) for f in (6.0, 5.0, 1.0)],
            [zeros(1, 1) for _ in 1:3], [fill(10.0, 1) for _ in 1:3],
            [zeros(1) for _ in 1:3], ["Y1"])
        for (T, t0) in ((5.0, 0.0), (4.0, 2.0))
            result = solve_coalescent_odes(tfgy, [4.0]; T_sample=T, t0, dt=0.3)
            @test all(isapprox(u[1], 4 / (1 + 4 / 100 * ((T + 1) * s - s^2 / 2));
                atol=1e-8) for (s, u) in zip(result.sol.t, result.sol.u))
            @test result.summary.t_end ≈ t0
        end
    end

    @testset "Three demes, group independence and relabeling" begin
        F = [3.0 1.0 0.5; 0.2 2.0 0.7; 0.4 0.1 1.0]
        G = [0.0 1.0 0.3; 0.2 0.0 0.1; 0.5 0.2 0.0]
        Y, samples = [50.0, 40.0, 30.0], [6.0, 4.0, 3.0]
        tfgy = constant_moment_tfgy(F, G, Y)
        groups = (:jumps, :coalescences, :branch_lengths, :cherries, :sackin,
            :fitch, :sankoff, :parsimony_lengths)
        C = [0 1 2; 3 0 1; 2 3 0]
        minimal = solve_coalescent_odes(tfgy, samples)
        full = solve_coalescent_odes(tfgy, samples; states=groups, costs=C)
        @test minimal.layout.nstates == 3
        @test length(full.layout.programs[:fitch].messages) == 7
        @test full.summary.lineages ≈ minimal.summary.lineages atol=1e-8
        @test sum(full.summary.descendants) ≈ sum(samples) atol=1e-8
        @test sum(full.summary.coalescences) ≈ sum(samples) - sum(full.summary.lineages) atol=1e-8
        for name in (:fitch, :sankoff)
            @test vec(sum(getproperty(full.summary, name); dims=2)) ≈ full.summary.lineages atol=1e-8
            @test sum(getproperty(full.summary, Symbol(name, :_lengths))) ≈
                sum(full.summary.branch_lengths) atol=1e-8
        end
        @test all(full.summary.external .>= -1e-10)
        @test all(vec(sum(full.summary.external; dims=2)) .<= full.summary.lineages .+ 1e-8)
        @test issymmetric(full.summary.cherries)
        @test full.summary.sankoff_score >= full.summary.fitch_score > 0

        permutation = [3, 1, 2]
        permuted = solve_coalescent_odes(
            constant_moment_tfgy(F[permutation, permutation], G[permutation, permutation], Y[permutation]),
            samples[permutation]; states=groups, character_states=permutation, costs=C)
        @test permuted.summary.lineages ≈ full.summary.lineages[permutation] atol=1e-8
        @test permuted.summary.jumps ≈ full.summary.jumps[permutation, permutation] atol=1e-8
        @test permuted.summary.cherries ≈ full.summary.cherries atol=1e-8
        @test permuted.summary.sankoff_score ≈ full.summary.sankoff_score atol=1e-8
        grouped = solve_coalescent_odes(tfgy, samples;
            states=(:fitch, :camin_sokal), character_states=[1, 2, 2])
        @test length(grouped.layout.programs[:fitch].messages) == 3
        @test grouped.summary.camin_sokal_12_score >= grouped.summary.fitch_score
    end

    @testset "State-budget and input failures" begin
        tfgy = constant_moment_tfgy(ones(3, 3), zeros(3, 3), fill(20.0, 3))
        @test_throws ArgumentError solve_coalescent_odes(tfgy, [2.0, 3.0, 4.0]; states=:unknown)
        @test_throws ArgumentError solve_coalescent_odes(tfgy, [2.0, 3.0, 4.0]; states=:parsimony_lengths)
        @test_throws ArgumentError solve_coalescent_odes(tfgy, [2.0, 3.0, 4.0];
            states=:fitch, max_parsimony_classes=6)
        @test_throws ArgumentError solve_coalescent_odes(tfgy, [2.0, 3.0, 4.0];
            states=:fitch, max_states=20)
        @test_throws ArgumentError solve_coalescent_odes(tfgy, [2.0, 3.0, 4.0]; states=:camin_sokal)
        @test_throws ArgumentError solve_coalescent_odes(tfgy, [2.0, 3.0, 4.0];
            states=:sankoff, costs=[0.0 0.5 1; 1 0 1; 1 1 0])
        @test_throws ArgumentError solve_coalescent_odes(tfgy, [2.0, 3.0, 4.0];
            states=:sankoff, costs=[0.0 Inf Inf; Inf 0 Inf; Inf Inf 0])
        @test_throws ArgumentError solve_coalescent_odes(tfgy, [21.0, 3.0, 4.0])
        @test_throws ArgumentError solve_coalescent_odes(tfgy, [0.0, 0.0, 0.0])
        @test_throws ArgumentError solve_coalescent_odes(tfgy, [2.0, 3.0])
        @test_throws ArgumentError solve_coalescent_odes(tfgy, [2.0, 3.0, 4.0]; dt=0)
        @test_throws ArgumentError solve_coalescent_odes(tfgy, [2.0, 3.0, 4.0];
            character_states=[1, 3, 3])
        zero_origin = TFGY([2.0, 0.0], tfgy.Fs, tfgy.Gs,
            [fill(20.0, 3), [20.0, 0.0, 20.0]], tfgy.deaths, tfgy.demes)
        @test_throws ArgumentError solve_coalescent_odes(zero_origin, [2.0, 3.0, 4.0])
        @test solve_coalescent_odes(zero_origin, [2.0, 3.0, 4.0]; t0=0.5).summary.t_end ≈ 0.5
        inconsistent = TFGY([2.0, 0.0], [zeros(1, 1), zeros(1, 1)],
            [zeros(1, 1), zeros(1, 1)], [[10.0], [1.0]], [zeros(1), zeros(1)], ["Y1"])
        @test_throws DomainError solve_coalescent_odes(inconsistent, [4.0])
    end

    @testset "ReactionSystem entry point" begin
        rn = @reaction_network begin
            birth, Z --> 2Z
        end
        result = solve_coalescent_odes(rn, [rn.Z], [rn.Z => 10.0], (0.0, 1.0),
            [rn.birth => 0.1], [rn.Z => 3]; states=(:branch_lengths, :sackin), res=100)
        @test result.summary.samples == [3.0]
        @test 1 < only(result.summary.lineages) < 3
        @test only(result.summary.descendants) ≈ 3.0
        @test_throws ArgumentError solve_coalescent_odes(rn, [rn.Z], [rn.Z => 10.0], (0.0, 1.0),
            [rn.birth => 0.1], [rn.Z => 2, rn.Z => 1])
    end
end

@testset "Compiled Sankoff automata versus ordinary tree dynamic programming" begin
    function ordinary_costs(tree, C)
        k = size(C, 1)
        tree isa Int && return [i == tree ? 0.0 : Inf for i in 1:k]
        a, b = ordinary_costs(tree[1], C), ordinary_costs(tree[2], C)
        return [minimum(C[i, j] + a[j] for j in 1:k) +
            minimum(C[i, j] + b[j] for j in 1:k) for i in 1:k]
    end
    function compiled_score(tree, program)
        tree isa Int && return (program.leaf_classes[tree], 0)
        a, score_a = compiled_score(tree[1], program)
        b, score_b = compiled_score(tree[2], program)
        return (program.parent[a, b], score_a + score_b + program.increments[a, b])
    end
    for C in ([0 1 1; 1 0 1; 1 1 0], [0 1 2; 3 0 1; 2 3 0],
              [0 3 1; 1 0 1; 1 1 0],
              [0.0 1.0 2.0; Inf 0.0 1.0; Inf Inf 0.0])
        program = Phylodynamics._compile_sankoff(C; max_classes=256)
        for a in 1:3, b in 1:3, c in 1:3, d in 1:3
            for tree in (((a, b), (c, d)), (((a, b), c), d))
                @test last(compiled_score(tree, program)) == minimum(ordinary_costs(tree, C))
            end
        end
    end
end

@testset "Corrected binary directional classes conserve true lineages" begin
    z = zeros(37)
    z[1:2] .= 1
    z[[15, 18, 19, 22, 23, 27]] .= 1
    du = similar(z)
    p = (beta=0.5, epsilon=1.0, phi=0.0, psi=0.0)
    Phylodynamics.Zoonosis._zoonosis_coalescent_rhs!(du, z, p, [0.0, 2.0, 1.0, 2.0])
    @test du[15] + du[16] ≈ du[1] ≈ 0.25
    @test du[17] + du[18] ≈ du[2] ≈ -0.5
    @test du[19] + du[20] ≈ du[1]
    @test du[21] + du[22] ≈ du[2]
    for (offset, weight) in ((15, 0.0), (19, 1.0))
        y = copy(z)
        y[offset:(offset + 3)] .= [1.0, 2.0, 3.0, 4.0]
        y[1:2] .= [3.0, 7.0]
        Phylodynamics.Zoonosis._zoonosis_coalescent_rhs!(du, y, p, [0.0, 2.0, 1.0, 2.0])
        derived = Phylodynamics.Zoonosis._general_parsimony_4state_derivs(
            1.0, 2.0, 3.0, 4.0, 2.0, 2.0, 0.0, 1.0, 0.0, 0.0, weight)
        @test collect(derived[1:4]) ≈ du[offset:(offset + 3)]
    end
    y = copy(z)
    y[23:28] .= [0.1, 0.3, 0.2, 0.4, 0.2, 0.1]
    y[1:2] .= [0.6, 0.7]
    Phylodynamics.Zoonosis._zoonosis_coalescent_rhs!(du, y, p, [0.0, 2.0, 1.0, 2.0])
    derived = Phylodynamics.Zoonosis._fitch_parsimony_6state_derivs(
        y[23:28]..., 2.0, 2.0, 0.0, 1.0, 0.0, 0.0)
    @test collect(derived[1:6]) ≈ du[23:28]
    @test derived[7] ≈ du[33]
    result = solve_coalescent_odes(model3_params(); n1=40.0, n2=10.0, dt=0.05)
    @test result.traj.A11 + result.traj.A12 ≈ result.traj.A1 atol=1e-7
    @test result.traj.A21 + result.traj.A22 ≈ result.traj.A2 atol=1e-7
    @test result.traj.B11 + result.traj.B12 ≈ result.traj.A1 atol=1e-7
    @test result.traj.B21 + result.traj.B22 ≈ result.traj.A2 atol=1e-7
end

@testset "General moment integration versus corrected two-host solver" begin
    for parameters in (model1_params(), model2_params(), model3_params())
        old = solve_coalescent_odes(parameters; n1=40.0, n2=10.0, dt=0.05)
        times = collect(range(12.85, 0.01; length=4000))
        sizes = [old.fwd.sol(t) for t in times]
        F = [[parameters.beta * u[1] * u[2] parameters.epsilon * parameters.beta * u[3] * u[2];
              parameters.psi * parameters.beta * u[1] * u[4] parameters.phi * parameters.beta * u[3] * u[4]]
             for u in sizes]
        tfgy = TFGY(times, F, [zeros(2, 2) for _ in times],
            [[u[2], u[4]] for u in sizes], [zeros(2) for _ in times], ["Y1", "Y2"])
        new = solve_coalescent_odes(tfgy, [40.0, 10.0];
            states=(:fitch, :camin_sokal), dt=0.05)
        # Compare trajectories before the continuous root stop, not the legacy
        # nearest-grid root summary.
        for row in eachrow(old.traj)
            row.s <= last(new.sol.t) || continue
            u = new.sol(row.s)
            @test u[new.layout.groups[:lineages]] ≈ [row.A1, row.A2] atol=2e-5
            for (name, expected) in ((:camin_sokal_12_score, row.J_pars_1to2),
                                     (:camin_sokal_21_score, row.J_pars_2to1),
                                     (:fitch_score, row.J_pars_fitch))
                @test only(u[new.layout.groups[name]]) ≈ expected atol=2e-5
            end
        end
    end
end
