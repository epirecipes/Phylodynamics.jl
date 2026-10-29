"""
    simulate_remaster_tree(rs::ReactionSystem, infected, u0, tspan, p, samples;
        remove_singletons=true, max_attempts=200, rng=Random.default_rng())

Simulate the supplied reaction system by direct SSA and trace sampled lineages
backwards through its transmission and physical-migration events. Supports any
number of infected demes, in the supplied species order, and homochronous samples
at `tspan[2]`. `samples` uses the same pair-map or `SampleConfiguration` interface
as [`TreeProblem`](@ref). Retry trajectories until all requested counts are
available; counts are never silently clamped.

Return `(tree, full_tree, binary_tree, ssa, jump_counts, jump_times)`.
The last two fields record forward-direction transitions for every deme pair.
`full_tree` retains singleton state changes; `binary_tree` collapses them.
If multiple sampled founder lineages remain at `tspan[1]`, fail explicitly rather
than fabricate their ancestry. The direct SSA requires time-homogeneous rates.
"""
function simulate_remaster_tree(rs::ReactionSystem, infected::AbstractVector,
    u0, tspan::Tuple{<:Real,<:Real}, p, samples;
    remove_singletons::Bool=true, max_attempts::Int=200,
    rng::AbstractRNG=Random.default_rng())
    inf_idx = _infected_indices(rs, infected)
    m = length(inf_idx)
    demes = _species_name.(infected)
    counts = _homochronous_sample_counts(samples, demes, tspan[2])
    ssa = simulate_tracked_ssa(rs, infected, u0, p, tspan;
        min_final_infected=counts, max_attempts=max_attempts, rng=rng)
    n_tips = sum(counts)
    tip_state = [k for k in 1:m for _ in 1:counts[k]]
    tip_label = ["tip_$(i)_$(demes[k])" for (i, k) in enumerate(tip_state)]
    node_state = copy(tip_state)
    node_time = fill(Float64(tspan[2]), n_tips)
    extant = [findall(==(k), tip_state) for k in 1:m]
    parents = Int[]
    children = Int[]
    lengths = Float64[]
    jump_counts = zeros(Int, m, m)
    jump_times = [Float64[] for _ in 1:m, _ in 1:m]

    function ancestor!(deme, time, descendants...)
        push!(node_state, deme)
        push!(node_time, time)
        node = length(node_state)
        for child in descendants
            push!(parents, node)
            push!(children, child)
            push!(lengths, max(1e-8, node_time[child] - time))
        end
        push!(extant[deme], node)
    end

    for event in Iterators.reverse(ssa.rx_events)
        sum(length, extant) <= 1 && break
        spec = ssa.specs[event.rx_idx]
        k, l = spec.parent_deme, spec.child_deme
        Ak, Al = length(extant[k]), length(extant[l])
        Yk, Yl = event.Y_after[k], event.Y_after[l]
        all(length(extant[i]) <= event.Y_after[i] for i in 1:m) ||
            throw(DomainError(event.Y_after, "ancestral counts exceed recorded populations"))
        if spec.role === :transmission && k == l
            if Ak >= 2 && rand(rng) < Ak * (Ak - 1.0) / (Yk * (Yk - 1.0))
                first = popat!(extant[k], rand(rng, 1:Ak))
                second = popat!(extant[k], rand(rng, 1:(Ak - 1)))
                ancestor!(k, event.time, first, second)
            end
        elseif Al > 0 && rand(rng) < Al / Yl
            child = popat!(extant[l], rand(rng, 1:Al))
            jump_counts[k, l] += 1
            push!(jump_times[k, l], event.time)
            if spec.role === :transmission && Ak > 0 && rand(rng) < Ak / Yk
                parent = popat!(extant[k], rand(rng, 1:Ak))
                ancestor!(k, event.time, parent, child)
            else
                ancestor!(k, event.time, child)
            end
        end
    end

    founders = reduce(vcat, extant)
    length(founders) == 1 || throw(ArgumentError(
        "$(length(founders)) sampled founder lineages remain at tspan[1]; " *
        "their common ancestry is not specified by the reaction history"))
    root = only(founders)
    internals = vcat(root, filter(!=(root), collect((n_tips + 1):length(node_state))))
    order = vcat(collect(1:n_tips), internals)
    index = invperm(order)
    jump_times .= sort.(jump_times)
    full_tree = PhyloTree(
        hcat(index[parents], index[children]), lengths, tip_label, tip_state,
        node_state[order], node_time[order], length(internals),
        m >= 2 ? jump_counts[1, 2] : 0, m >= 2 ? jump_counts[2, 1] : 0,
        m >= 2 ? copy(jump_times[1, 2]) : Float64[],
        m >= 2 ? copy(jump_times[2, 1]) : Float64[])
    binary_tree = collapse_singletons(full_tree)
    return (tree=remove_singletons ? binary_tree : full_tree,
        full_tree=full_tree, binary_tree=binary_tree, ssa=ssa,
        jump_counts=jump_counts, jump_times=jump_times)
end

"""Alias for [`simulate_remaster_tree`](@ref); accepts the same arguments and keywords."""
const simulate_forward_backward_tree = simulate_remaster_tree
