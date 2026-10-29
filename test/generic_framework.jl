module GenericFrameworkTests
using Test
using Phylodynamics
using Catalyst
using ModelingToolkit
using SciMLBase
using StableRNGs

@testset "Generic framework and explicit study namespace" begin
    @test parentmodule(Phylodynamics.Zoonosis.ZoonosisParams) === Phylodynamics.Zoonosis
    for name in (:ZoonosisParams, :model1_params, :build_zoonosis_reaction_system,
        :compute_canonical_summary, :sweep_psi, :generate_all_figures,
        :fit_structured_coalescent, :zoonosis_demographic_process)
        @test !isdefined(Phylodynamics, name)
        @test name in names(Phylodynamics.Zoonosis)
    end

    t = ModelingToolkit.t_nounits
    @parameters growth transfer
    @species A(t) B(t) C(t)
    @named one = ReactionSystem([Reaction(growth, [A], [A], [1], [2])],
        t, [A], [growth]; combinatoric_ratelaws=false)
    one = complete(one)
    single = simulate_remaster_tree(one, [A], [A => 1], (0.0, 4.0),
        [growth => 1.0], [A => 4]; rng=StableRNG(101))
    @test ntips(single.tree) == 4
    @test single.tree.nnode == 3
    @test all(==(1), single.full_tree.node_state)
    @test single.jump_counts == zeros(Int, 1, 1)
    @test all(length(event.Y_after) == 1 for event in single.ssa.rx_events)
    @test single.ssa.times[end] == 4.0
    @test_throws ArgumentError simulate_remaster_tree(one, [A], [A => 3],
        (0.0, 1.0), [growth => 0.0], [A => 3])
    @test_throws ArgumentError simulate_forward_master_tree(one, [A], [A => 3],
        (0.0, 1.0), [growth => 0.0], [A => 3])
    @test_throws ArgumentError simulate_tracked_ssa(one, [A], [A => 1.5],
        [growth => 1.0], (0.0, 1.0))
    @test_throws ArgumentError simulate_tracked_ssa(one, [A], [A => 1],
        [growth => 1.0], (0.0, 1.0); min_final_infected=[1, 0])
    @test_throws ArgumentError simulate_remaster_tree(one, [A], [A => 1],
        (0.0, 1.0), [growth => 1.0], [B => 2])
    serial = SampleConfiguration(["A", "A"], [0.5, 1.0])
    @test_throws ArgumentError simulate_remaster_tree(one, [A], [A => 1],
        (0.0, 1.0), [growth => 1.0], serial)

    @named chain = ReactionSystem([
        Reaction(growth, [A], [A], [1], [2]),
        Reaction(growth, [B], [B], [1], [2]),
        Reaction(growth, [C], [C], [1], [2]),
        Reaction(transfer, [A], [B], [1], [1]),
        Reaction(transfer, [B], [C], [1], [1])],
        t, [A, B, C], [growth, transfer]; combinatoric_ratelaws=false)
    chain = complete(chain)
    infected = [B, C, A]
    u0 = [A => 1, B => 0, C => 0]
    p = [growth => 1.5, transfer => 0.7]
    samples = [A => 3, B => 3, C => 3]
    result = simulate_remaster_tree(chain, infected, u0, (0.0, 3.0),
        p, samples; rng=StableRNG(77))
    @test ntips(result.tree) == 9
    @test result.tree.nnode == 8
    @test [count(==(k), result.tree.tip_state) for k in 1:3] == [3, 3, 3]
    @test result.full_tree.node_state[root_node(result.full_tree)] == 3
    @test result.jump_counts[3, 1] > 0
    @test result.jump_counts[1, 2] > 0
    @test sum(result.jump_counts) == result.jump_counts[3, 1] + result.jump_counts[1, 2]
    @test length.(result.jump_times) == result.jump_counts
    @test any(result.ssa.specs[e.rx_idx].role === :migration for e in result.ssa.rx_events)
    @test all(length(e.Y_after) == 3 for e in result.ssa.rx_events)
    @test all(x -> x in 1:3, result.full_tree.node_state)
    @test all(label -> occursin(r"_(A|B|C)$", label), result.tree.tip_label)
    restored = read_newick(write_newick(result.full_tree; include_metadata=true))
    @test Dict(restored.tip_label .=> restored.tip_state) ==
        Dict(result.full_tree.tip_label .=> result.full_tree.tip_state)
    @test sort(restored.node_state) == sort(result.full_tree.node_state)

    problem = TreeProblem(chain, infected, u0, (0.0, 3.0), p, samples)
    tree = solve(problem, ReMASTERStepper(); rng=StableRNG(77))
    @test tree.edge == result.tree.edge
    @test tree.tip_state == result.tree.tip_state
    @test tree.node_time == result.tree.node_time

    model = ModelFGY(demes=["deme_b", "deme_a", "deme_c"], tfin=2.0)
    events = [Event(SAMPLE, 0.0, "deme_a", "deme_a"),
        Event(SAMPLE, 0.0, "deme_b", "deme_b"),
        Event(MIGRATION, 0.5, "deme_b", "deme_a"),
        Event(COALESCENT, 1.0, "deme_b", "deme_b")]
    simulated = SimTree(events, model; record_events=true, rng=StableRNG(11))
    @test simulated.event_history !== nothing
    @test [event.height for event in simulated.event_history] == [0.0, 0.0, 0.5, 1.0]
    converted = simtree_to_phylotree(simulated; T_sample=2.0)
    @test converted.tip_state == [2, 1]
    @test converted.true_jumps_12 == 1
    @test converted.true_jumps_21 == 0
    @test !any(endswith(label, "_I") || endswith(label, "_Y") for label in converted.tip_label)

    forest_events = [Event(SAMPLE, 0.0, "deme_b", "deme_b"),
        Event(SAMPLE, 0.0, "deme_b", "deme_b"),
        Event(SAMPLE, 0.0, "deme_a", "deme_a"),
        Event(COALESCENT, 0.5, "deme_b", "deme_b")]
    forest = SimTree(forest_events, model; record_events=true, rng=StableRNG(12))
    @test count(e -> e.type == COALESCENT, forest.event_history) == 1
    @test forest.nNode == 2
    @test maximum(e.height for e in forest.event_history) == 0.5
    @test maximum(Float64.(forest.heights)) > 0.5

    @named varying = ReactionSystem([Reaction(growth * t, [A], [A], [1], [2])],
        t, [A], [growth]; combinatoric_ratelaws=false)
    varying = complete(varying)
    @test_throws ArgumentError simulate_tracked_ssa(varying, [A], [A => 1],
        [growth => 1.0], (0.0, 1.0))
end
end
