# ==============================================================================
# reaction_network.jl
# Requirement (a) & (b):
# - Defines epidemiological reaction networks using ModelingToolkit.jl
#   (@independent_variables, @parameters, @species, Reaction, ReactionSystem)
#   and Catalyst.jl (@reaction_network)
# - Automatically classifies reactions into :transmission (F_{kl}), :migration (G_{kl}),
#   :removal (D_k), and :demographic directly from stoichiometry on infected species
# - Converts ReactionSystem to:
#   * ODESystem / ODEProblem (solved with OrdinaryDiffEqTsit5 / Vern7 / Rosenbrock23)
#   * JumpSystem / JumpProblem (solved with JumpProcesses.jl Direct() + SSAStepper())
#   * Coalescent.jl ModelFGY and phydynR TFGY
# - Implements ReMASTER-style forward-backward stochastic reaction tree generation
#   (producing both full singleton-annotated reaction trees and collapsed binary trees)
#   and deterministic Coalescent.jl structured coalescent tree simulation.
# ==============================================================================

"""
    PhyloReactionSpec

Classification of a single `ModelingToolkit.Reaction` in an epidemiological `ReactionSystem`
with respect to the vector of infected species `infected_species = [Y1, Y2, ...]`.

- `role`: `:transmission` (`parent_deme -> child_deme`, net infected `+1`),
          `:migration` (`parent_deme -> child_deme`, net infected `0`, `k != l`),
          `:removal` (`parent_deme`, net infected `-1`), or
          `:demographic` (uninfected compartments only).
"""
struct PhyloReactionSpec
    rx_index::Int
    role::Symbol
    parent_deme::Int
    child_deme::Int
    rx::Catalyst.Reaction
end

"""
    classify_phylodynamic_reactions(rs::ReactionSystem, infected_species::AbstractVector) -> Vector{PhyloReactionSpec}

Inspect the stoichiometry of each `Reaction` in `rs` on `infected_species` and classify it
into `:transmission` (`k -> l`), `:migration` (`k -> l`), `:removal` (`k -> 0`), or `:demographic`.

Just like ReMASTER and `ReproductiveNumbers.jl`, a reaction is:
- `:transmission` from deme `k` to deme `l` if `k` is an infected substrate (stoich 1)
  and the infected products are `(k, l)` (net infected stoichiometry `+1`).
- `:migration` from deme `k` to deme `l` (`k != l`) if `k` is an infected substrate (stoich 1)
  and `l` is the single infected product (net infected stoichiometry `0`).
- `:removal` from deme `k` if `k` is an infected substrate (stoich 1) and there are no
  infected products (net infected stoichiometry `-1`).
"""
function classify_phylodynamic_reactions(rs::ReactionSystem, infected_species::AbstractVector)
    rxs = Catalyst.reactions(rs)
    inf_syms = Num.(infected_species)
    specs = PhyloReactionSpec[]
    sizehint!(specs, length(rxs))

    for (idx, rx) in enumerate(rxs)
        sub_demes = Int[]
        for (s, st) in zip(rx.substrates, rx.substoich)
            k = findfirst(x -> isequal(Num(s), x), inf_syms)
            if k !== nothing
                for _ in 1:Int(st)
                    push!(sub_demes, k)
                end
            end
        end
        prod_demes = Int[]
        for (p, pt) in zip(rx.products, rx.prodstoich)
            k = findfirst(x -> isequal(Num(p), x), inf_syms)
            if k !== nothing
                for _ in 1:Int(pt)
                    push!(prod_demes, k)
                end
            end
        end

        if length(sub_demes) == 1 && length(prod_demes) == 2
            parent_k = sub_demes[1]
            # The newly produced infected individual is the remaining product after removing one parent_k
            if prod_demes[1] == parent_k
                child_l = prod_demes[2]
            elseif prod_demes[2] == parent_k
                child_l = prod_demes[1]
            else
                child_l = prod_demes[1]
            end
            push!(specs, PhyloReactionSpec(idx, :transmission, parent_k, child_l, rx))
        elseif length(sub_demes) == 1 && length(prod_demes) == 1 && sub_demes[1] != prod_demes[1]
            push!(specs, PhyloReactionSpec(idx, :migration, sub_demes[1], prod_demes[1], rx))
        elseif length(sub_demes) == 1 && isempty(prod_demes)
            push!(specs, PhyloReactionSpec(idx, :removal, sub_demes[1], 0, rx))
        else
            push!(specs, PhyloReactionSpec(idx, :demographic, 0, 0, rx))
        end
    end
    return specs
end

"""
    reaction_system_to_fgy(rs::ReactionSystem, infected_species::AbstractVector, u0_map, p_map, tspan; res::Int=1000, solver=Tsit5()) -> TFGY

Convert any `ModelingToolkit.ReactionSystem` and its classified reactions into a `phydynR` `TFGY`
trajectory `(times, Fs, Gs, Ys, deaths, demes)` by solving the MTK `ODESystem` forward and
evaluating each reaction's `oderatelaw(rx)` along the reverse-time grid.
"""
function reaction_system_to_fgy(
    rs::ReactionSystem,
    infected_species::AbstractVector,
    u0_map,
    p_map,
    tspan::Tuple{Float64,Float64};
    res::Int = 1000,
    solver = Tsit5()
)::TFGY
    specs = classify_phylodynamic_reactions(rs, infected_species)
    odesys = complete(Catalyst.ode_model(rs))
    saveat = collect(range(tspan[1], tspan[2]; length=res))
    oprob = ODEProblem(odesys, u0_map, tspan, p_map)
    sol = solve(oprob, solver; saveat=saveat, reltol=1e-9, abstol=1e-11)
    SciMLBase.successful_retcode(sol) || error("forward demographic integration failed: $(sol.retcode)")

    m = length(infected_species)
    demes = [String(Symbolics.tosymbol(s; escape=false)) for s in infected_species]
    comb = Catalyst.get_combinatoric_ratelaws(rs)
    rate_laws = [Num(Catalyst.oderatelaw(sp.rx; combinatoric_ratelaw=comb)) for sp in specs]
    spcs = Catalyst.species(rs)
    ps = Catalyst.parameters(rs)
    t_var = ModelingToolkit.get_iv(rs)
    compiled_rates = [Symbolics.build_function(rl, spcs, ps, t_var; expression=Val{false}) for rl in rate_laws]

    p_dict = Dict(Num(k) => Float64(v) for (k, v) in p_map)
    p_vec = [p_dict[Num(pr)] for pr in ps]
    sp_vecs = [sol[s] for s in spcs]
    inf_vecs = [sol[s] for s in infected_species]

    times_rev = reverse(sol.t)
    nt = length(times_rev)
    Fs = Vector{Matrix{Float64}}(undef, nt)
    Gs = Vector{Matrix{Float64}}(undef, nt)
    Ys = Vector{Vector{Float64}}(undef, nt)
    deaths = Vector{Vector{Float64}}(undef, nt)
    u_buf = zeros(Float64, length(spcs))

    for idx in 1:nt
        fwd_i = nt - idx + 1
        t_val = sol.t[fwd_i]
        for j in eachindex(spcs)
            u_buf[j] = max(0.0, sp_vecs[j][fwd_i])
        end
        F = zeros(Float64, m, m)
        G = zeros(Float64, m, m)
        Y = [max(0.0, inf_vecs[k][fwd_i]) for k in 1:m]
        D = zeros(Float64, m)
        for (r_i, sp) in enumerate(specs)
            sp.role === :demographic && continue
            r_val = max(0.0, Float64(compiled_rates[r_i](u_buf, p_vec, t_val)))
            if sp.role === :transmission
                F[sp.parent_deme, sp.child_deme] += r_val
            elseif sp.role === :migration
                G[sp.parent_deme, sp.child_deme] += r_val
            elseif sp.role === :removal
                D[sp.parent_deme] += r_val
            end
        end
        Fs[idx] = F
        Gs[idx] = G
        Ys[idx] = Y
        deaths[idx] = D
    end

    return TFGY(times_rev, Fs, Gs, Ys, deaths, demes)
end

# Event record for ReMASTER-style forward-backward stochastic reaction tree simulation
struct _RxEventRecord
    time::Float64
    rx_idx::Int
    Y_after::Vector{Int}
end

function _infected_indices(rs, infected_species)
    isempty(infected_species) && throw(ArgumentError("at least one infected species is required"))
    indices = Int[]
    spcs = Catalyst.species(rs)
    for infected in infected_species
        index = findfirst(s -> isequal(Num(s), Num(infected)), spcs)
        index === nothing && throw(ArgumentError("infected species $infected is not in the reaction system"))
        index in indices && throw(ArgumentError("infected species must be distinct"))
        push!(indices, index)
    end
    return indices
end

function _homochronous_sample_counts(samples, demes, t_end)
    scfg = samples isa SampleConfiguration ? samples : SampleConfiguration(samples, t_end)
    for (deme, time) in scfg.sconf
        deme in demes || throw(ArgumentError("sample deme $deme is not in the infected species"))
        time == t_end || throw(ArgumentError("SSA tree sampling must be contemporaneous at tspan[2]"))
    end
    counts = [count(x -> x[1] == deme, scfg.sconf) for deme in demes]
    sum(counts) >= 2 || throw(ArgumentError("at least two sampled tips are required"))
    return counts
end

"""
    simulate_tracked_ssa(rs::ReactionSystem, infected_species, u0_map, p_map, tspan; min_final_infected=nothing, max_attempts=200, rng=Random.default_rng())

Run an exact Gillespie SSA trajectory directly from the `ModelingToolkit.ReactionSystem`'s
`jumpratelaw`s and `netstoichmat(rs)`, recording every transmission and migration reaction
event `(time, rx_idx, Y_after)` needed for ReMASTER's backward-in-time tree construction.
`Y_after` contains all infected populations in the supplied species order.
Initial counts must be nonnegative integers. Optional `min_final_infected` supplies
one minimum per deme; the default imposes no survival conditioning.
Explicitly time-dependent propensities are not supported by this direct SSA.
"""
function simulate_tracked_ssa(
    rs::ReactionSystem,
    infected_species::AbstractVector,
    u0_map,
    p_map,
    tspan::Tuple{<:Real,<:Real};
    min_final_infected = nothing,
    max_attempts::Int = 200,
    rng::AbstractRNG = Random.default_rng()
)
    inf_idx = _infected_indices(rs, infected_species)
    m = length(inf_idx)
    min_counts = min_final_infected === nothing ? zeros(Int, m) : collect(min_final_infected)
    length(min_counts) == m && all(x -> x isa Integer && x >= 0, min_counts) ||
        throw(ArgumentError("min_final_infected must contain one nonnegative integer per deme"))
    max_attempts > 0 || throw(ArgumentError("max_attempts must be positive"))
    all(isfinite, tspan) && tspan[1] < tspan[2] ||
        throw(ArgumentError("tspan must be finite and increasing"))
    all(pair -> isfinite(last(pair)) && last(pair) >= 0 && isinteger(last(pair)), u0_map) ||
        throw(ArgumentError("initial population counts must be nonnegative integers"))
    specs = classify_phylodynamic_reactions(rs, infected_species)
    spcs = Catalyst.species(rs)
    ps = Catalyst.parameters(rs)
    t_var = ModelingToolkit.get_iv(rs)
    comb = Catalyst.get_combinatoric_ratelaws(rs)
    jump_laws = [Num(Catalyst.jumpratelaw(sp.rx; combinatoric_ratelaw=comb)) for sp in specs]
    any(law -> any(v -> isequal(v, t_var), Symbolics.get_variables(law)), jump_laws) &&
        throw(ArgumentError("time-dependent propensities are not supported by the direct SSA"))
    compiled_props = [Symbolics.build_function(jl, spcs, ps, t_var; expression=Val{false}) for jl in jump_laws]
    net_stoich = Catalyst.netstoichmat(rs)

    u0_dict = Dict(Num(k) => round(Int, v) for (k, v) in u0_map)
    p_dict = Dict(Num(k) => Float64(v) for (k, v) in p_map)
    u0_vec = [u0_dict[Num(s)] for s in spcs]
    p_vec = [p_dict[Num(pr)] for pr in ps]
    nrx = length(specs)

    props = zeros(Float64, nrx)
    u = similar(u0_vec)

    for attempt in 1:max_attempts
        copyto!(u, u0_vec)
        t = Float64(tspan[1])
        t_end = Float64(tspan[2])
        t_hist = Float64[t]
        u_hist = Vector{Int}[copy(u)]
        rx_events = _RxEventRecord[]
        sizehint!(rx_events, 30000)

        while t < t_end
            tot_rate = 0.0
            @inbounds for r in 1:nrx
                val = Float64(compiled_props[r](u, p_vec, t))
                isfinite(val) && val >= 0 ||
                    throw(DomainError(val, "reaction $r has an invalid SSA propensity"))
                props[r] = val
                tot_rate += val
            end
            isfinite(tot_rate) || throw(DomainError(tot_rate, "total SSA propensity overflowed"))
            tot_rate <= 0.0 && break
            dt = randexp(rng) / tot_rate
            if t + dt > t_end
                t = t_end
                push!(t_hist, t)
                push!(u_hist, copy(u))
                break
            end
            t += dt
            # Select reaction proportional to propensity
            u_draw = rand(rng) * tot_rate
            cum = 0.0
            chosen = nrx
            @inbounds for r in 1:nrx
                cum += props[r]
                if u_draw <= cum
                    chosen = r
                    break
                end
            end
            @inbounds for s_i in eachindex(u)
                u[s_i] += net_stoich[s_i, chosen]
            end
            all(>=(0), u) || throw(DomainError(copy(u), "SSA reaction produced a negative population"))
            role = specs[chosen].role
            if role === :transmission || role === :migration
                push!(rx_events, _RxEventRecord(t, chosen, u[inf_idx]))
            end
            if length(t_hist) == 1 || (t - t_hist[end]) >= 0.05
                push!(t_hist, t)
                push!(u_hist, copy(u))
            end
        end

        if t_hist[end] < t_end
            push!(t_hist, t_end)
            push!(u_hist, copy(u))
        end
        if all(u[inf_idx] .>= min_counts)
            return (
                times = t_hist,
                states = u_hist,
                final_state = copy(u),
                rx_events = rx_events,
                specs = specs,
                attempt = attempt
            )
        end
    end
    error("simulate_tracked_ssa did not reach min_final_infected=$min_final_infected in $max_attempts attempts")
end

"""
    ModelFGY(rs::ReactionSystem, infected_species::AbstractVector, u0_map, tspan::Tuple{<:Real,<:Real}, p_map; modelname=String(nameof(rs))) -> ModelFGY

Construct a `Coalescent.jl` `ModelFGY` directly from a `ModelingToolkit.ReactionSystem`
(or `Catalyst.@reaction_network`), without any YAML or string parsing:
- Infected species `infected_species = [Y1, Y2, ...]` become the structured coalescent `demes`.
- Remaining species in `species(rs)` become `nondemes` (e.g., `X1, X2`).
- Each reaction's stoichiometry on `infected_species` classifies it into birth (`RXN_BIRTH`),
  migration (`RXN_MIG`), death (`RXN_DEATH`), or non-deme dynamics (`RXN_DYNVAR`), with
  compiled rate functions built via `Catalyst.oderatelaw` and `Symbolics.build_function`.
"""
function ModelFGY(
    rs::ReactionSystem,
    infected_species::AbstractVector,
    u0_map,
    tspan::Tuple{<:Real,<:Real},
    p_map;
    modelname::AbstractString = String(nameof(rs))
)::ModelFGY
    specs = classify_phylodynamic_reactions(rs, infected_species)
    spcs = Catalyst.species(rs)
    ps = Catalyst.parameters(rs)
    t_var = ModelingToolkit.get_iv(rs)
    comb = Catalyst.get_combinatoric_ratelaws(rs)

    demes = [_species_name(s) for s in infected_species]
    all_sp_names = [_species_name(s) for s in spcs]
    nondemes = [nm for nm in all_sp_names if !(nm in demes)]
    ordered_names = vcat(demes, nondemes)

    # Map ordered state vector u (ordered [demes..., nondemes...]) to Catalyst species order
    sp_pos_in_u = [findfirst(==(nm), ordered_names) for nm in all_sp_names]
    p_name_syms = [Symbol(_species_name(pr)) for pr in ps]

    birthrxn = CoalescentReaction[]
    migrationrxn = CoalescentReaction[]
    deathrxn = CoalescentReaction[]

    for sp in specs
        sp.role === :demographic && continue
        rl = Num(Catalyst.oderatelaw(sp.rx; combinatoric_ratelaw=comb))
        f_raw = Symbolics.build_function(rl, spcs, ps, t_var; expression=Val{false})
        rate_closure = let f = f_raw, pos = sp_pos_in_u, psyms = p_name_syms
            function (u, par_nt, t)
                u_sp = [max(0.0, Float64(u[idx])) for idx in pos]
                p_sp = [Float64(getproperty(par_nt, psym)) for psym in psyms]
                return max(0.0, Float64(f(u_sp, p_sp, t)))
            end
        end
        if sp.role === :transmission
            src = demes[sp.parent_deme]
            rec = demes[sp.child_deme]
            push!(birthrxn, CoalescentReaction(src, rec, RXN_BIRTH, :mtk_rate; rate_fn=rate_closure))
        elseif sp.role === :migration
            src = demes[sp.parent_deme]
            rec = demes[sp.child_deme]
            push!(migrationrxn, CoalescentReaction(src, rec, RXN_MIG, :mtk_rate; rate_fn=rate_closure))
        elseif sp.role === :removal
            src = demes[sp.parent_deme]
            push!(deathrxn, CoalescentReaction(src, RXN_DEATH, :mtk_rate; rate_fn=rate_closure))
        end
    end

    # Build exact symbolic net ODE right-hand sides for each non-deme variable from netstoichmat(rs)
    net_stoich = Catalyst.netstoichmat(rs)
    rxs = Catalyst.reactions(rs)
    rate_laws = [Num(Catalyst.oderatelaw(rx; combinatoric_ratelaw=comb)) for rx in rxs]
    nondemerxn = CoalescentReaction[]
    for nd_name in nondemes
        sp_idx = findfirst(==(nd_name), all_sp_names)
        rhs_sym = Num(0)
        for (k, rl) in enumerate(rate_laws)
            c = net_stoich[sp_idx, k]
            !iszero(c) && (rhs_sym += c * rl)
        end
        f_nd = Symbolics.build_function(rhs_sym, spcs, ps, t_var; expression=Val{false})
        nd_closure = let f = f_nd, pos = sp_pos_in_u, psyms = p_name_syms
            function (u, par_nt, t)
                u_sp = [max(0.0, Float64(u[idx])) for idx in pos]
                p_sp = [Float64(getproperty(par_nt, psym)) for psym in psyms]
                return Float64(f(u_sp, p_sp, t))
            end
        end
        push!(nondemerxn, CoalescentReaction(nd_name, RXN_DYNVAR, :mtk_ode; rate_fn=nd_closure))
    end

    initial = Dict(_species_name(k) => Float64(v) for (k, v) in u0_map)
    parameters = Dict(_species_name(k) => Float64(v) for (k, v) in p_map)

    return ModelFGY(;
        modelname = String(modelname),
        birthrxn = birthrxn,
        migrationrxn = migrationrxn,
        deathrxn = deathrxn,
        nondemerxn = nondemerxn,
        demes = demes,
        nondemes = nondemes,
        initial = initial,
        t0 = Float64(tspan[1]),
        tfin = Float64(tspan[2]),
        parameters = parameters
    )
end

"""
    build_demographic_process(rs::ReactionSystem, infected_species::AbstractVector, u0_map, tspan::Tuple{<:Real,<:Real}, p_map) -> Function

Build a `phydynR`-compatible demographic process callable `dm(theta, x0, t0, t1; res=1000)`
directly from a `ModelingToolkit.ReactionSystem`!
"""
function build_demographic_process(
    rs::ReactionSystem,
    infected_species::AbstractVector,
    u0_map,
    tspan::Tuple{<:Real,<:Real},
    p_map
)
    base_mod = ModelFGY(rs, infected_species, u0_map, tspan, p_map)
    return build_demographic_process(base_mod)
end

"""
    simulate_forward_master_tree(rs::ReactionSystem, infected_species::AbstractVector, u0_map, tspan::Tuple{<:Real,<:Real}, p_map, samples; remove_singletons::Bool=true, min_final_infected=nothing, max_attempts::Int=200, rng=Random.default_rng())

Simulate a genealogy **forward in time** directly from a `ModelingToolkit.ReactionSystem`
(like MASTER / Vaughan & Drummond 2013):
1. Tracks every individual infected host in deme `k` as a live node in `active_nodes[k]`
   during the forward Gillespie SSA loop (`length(active_nodes[k]) == Y_k(t)`).
2. When a transmission reaction `k -> l` fires at time `t`:
   - Picks a random parent lineage `u` uniformly from `active_nodes[k]`.
   - Creates an internal branching node `b` at time `t` in deme `k` (attached to `u`),
     spawning two daughter lineages: one continuing in `active_nodes[k]` and the newly
     infected offspring entering `active_nodes[l]` (via a state-transition node if `l != k`).
3. When a migration reaction `k -> l` (`k != l`) fires at time `t`:
   - Picks a random lineage `u` from `active_nodes[k]`, records a singleton migration node
     at time `t`, and moves it to `active_nodes[l]`.
4. When a removal/recovery/death reaction `k -> 0` fires at time `t`:
   - Picks a random lineage `u` from `active_nodes[k]` and terminates it (unsampled/extinct lineage).
5. At `T_sample = tspan[2]`, samples the requested number of lineages without replacement
   from `active_nodes[k]` and prunes all unsampled/extinct side branches to yield the exact
   sampled phylogeny (`full_tree` with singleton transition nodes and `binary_tree` with
   singletons collapsed).
"""
function simulate_forward_master_tree(
    rs::ReactionSystem,
    infected_species::AbstractVector,
    u0_map,
    tspan::Tuple{<:Real,<:Real},
    p_map,
    samples;
    remove_singletons::Bool = true,
    min_final_infected::Union{Nothing,Vector{Int}} = nothing,
    max_attempts::Int = 200,
    rng::AbstractRNG = Random.default_rng()
)
    specs = classify_phylodynamic_reactions(rs, infected_species)
    spcs = Catalyst.species(rs)
    ps = Catalyst.parameters(rs)
    t_var = ModelingToolkit.get_iv(rs)
    comb = Catalyst.get_combinatoric_ratelaws(rs)
    jump_laws = [Num(Catalyst.jumpratelaw(sp.rx; combinatoric_ratelaw=comb)) for sp in specs]
    compiled_props = [Symbolics.build_function(jl, spcs, ps, t_var; expression=Val{false}) for jl in jump_laws]
    net_stoich = Catalyst.netstoichmat(rs)

    u0_dict = Dict(Num(k) => round(Int, v) for (k, v) in u0_map)
    p_dict = Dict(Num(k) => Float64(v) for (k, v) in p_map)
    u0_vec = [u0_dict[Num(s)] for s in spcs]
    p_vec = [p_dict[Num(pr)] for pr in ps]
    m = length(infected_species)
    inf_idx = _infected_indices(rs, infected_species)
    demes = [_species_name(s) for s in infected_species]

    req_counts = _homochronous_sample_counts(samples, demes, tspan[2])
    min_req = min_final_infected === nothing ? req_counts : min_final_infected

    nrx = length(specs)
    props = zeros(Float64, nrx)
    u = similar(u0_vec)
    t_start = Float64(tspan[1])
    t_end = Float64(tspan[2])

    for attempt in 1:max_attempts
        copyto!(u, u0_vec)
        t = t_start

        # Forward phylogenetic graph storage:
        # Each created node has birth_time, deme_state, parent_node, and children
        f_parent = Int[]
        f_deme = Int[]
        f_time = Float64[]
        active_nodes = [Int[] for _ in 1:m]

        for k in 1:m
            init_cnt = u[inf_idx[k]]
            for _ in 1:init_cnt
                push!(f_parent, 0)
                push!(f_deme, k)
                push!(f_time, t_start)
                push!(active_nodes[k], length(f_parent))
            end
        end

        while t < t_end
            tot_rate = 0.0
            @inbounds for r in 1:nrx
                val = max(0.0, Float64(compiled_props[r](u, p_vec, t)))
                props[r] = val
                tot_rate += val
            end
            tot_rate <= 0.0 && break
            dt = randexp(rng) / tot_rate
            if t + dt > t_end
                t = t_end
                break
            end
            t += dt

            u_draw = rand(rng) * tot_rate
            cum = 0.0
            chosen = nrx
            @inbounds for r in 1:nrx
                cum += props[r]
                if u_draw <= cum
                    chosen = r
                    break
                end
            end

            @inbounds for s_i in eachindex(u)
                u[s_i] += net_stoich[s_i, chosen]
            end

            sp = specs[chosen]
            if sp.role === :transmission
                k = sp.parent_deme
                l = sp.child_deme
                pool_k = active_nodes[k]
                if !isempty(pool_k)
                    idx_p = rand(rng, 1:length(pool_k))
                    parent_node = pool_k[idx_p]
                    # Internal branching node at transmission time t in deme k
                    push!(f_parent, parent_node)
                    push!(f_deme, k)
                    push!(f_time, t)
                    branch_node = length(f_parent)
                    # Parent continues as daughter 1 in deme k
                    pool_k[idx_p] = branch_node
                    # Offspring enters deme l (if l != k, record state transition at time t)
                    push!(f_parent, branch_node)
                    push!(f_deme, l)
                    push!(f_time, t)
                    child_node = length(f_parent)
                    push!(active_nodes[l], child_node)
                end
            elseif sp.role === :migration
                k = sp.parent_deme
                l = sp.child_deme
                pool_k = active_nodes[k]
                if !isempty(pool_k)
                    idx_m = rand(rng, 1:length(pool_k))
                    old_node = pool_k[idx_m]
                    deleteat!(pool_k, idx_m)
                    push!(f_parent, old_node)
                    push!(f_deme, l)
                    push!(f_time, t)
                    mig_node = length(f_parent)
                    push!(active_nodes[l], mig_node)
                end
            elseif sp.role === :removal
                k = sp.parent_deme
                pool_k = active_nodes[k]
                if !isempty(pool_k)
                    idx_d = rand(rng, 1:length(pool_k))
                    deleteat!(pool_k, idx_d)
                end
            end
        end

        ok = all(k -> length(active_nodes[k]) >= min_req[k] || req_counts[k] == 0, 1:m)
        !ok && continue

        # Sample tips without replacement at t_end
        sampled_leaf_nodes = Int[]
        tip_label = String[]
        tip_state = Int[]
        for k in 1:m
            n_k = min(req_counts[k], length(active_nodes[k]))
            chosen_k = sample(rng, active_nodes[k], n_k; replace=false)
            suffix = k == 2 ? "Y" : "I"
            for (j, anc_node) in enumerate(chosen_k)
                push!(f_parent, anc_node)
                push!(f_deme, k)
                push!(f_time, t_end)
                leaf_id = length(f_parent)
                push!(sampled_leaf_nodes, leaf_id)
                push!(tip_label, "tip_$(suffix)_$(j)_$(suffix)")
                push!(tip_state, k)
            end
        end

        n_tips = length(sampled_leaf_nodes)
        n_tips < 2 && continue

        # Prune forward transmission graph to ancestors of sampled_leaf_nodes
        n_f = length(f_parent)
        sample_desc_count = zeros(Int, n_f)
        is_sampled_leaf = falses(n_f)
        for lf in sampled_leaf_nodes
            is_sampled_leaf[lf] = true
            curr = lf
            while curr > 0
                sample_desc_count[curr] += 1
                curr = f_parent[curr]
            end
        end

        # Build induced children lists on nodes with sample_desc_count > 0
        f_children = [Int[] for _ in 1:n_f]
        roots_used = Int[]
        for u_id in 1:n_f
            sample_desc_count[u_id] == 0 && continue
            p_id = f_parent[u_id]
            if p_id > 0 && sample_desc_count[p_id] > 0
                push!(f_children[p_id], u_id)
            else
                push!(roots_used, u_id)
            end
        end

        length(roots_used) == 1 || throw(ArgumentError(
            "$(length(roots_used)) sampled founder lineages remain at tspan[1]; " *
            "their common ancestry is not specified by the reaction history"))
        orig_root = roots_used[1]

        # Skip redundant degree-1 chain above MRCA unless state changed
        while !is_sampled_leaf[orig_root] && length(f_children[orig_root]) == 1
            orig_root = f_children[orig_root][1]
        end

        # Keep: (a) sampled leaves, (b) bifurcating internal nodes (length(f_children) >= 2),
        # and (c) singleton state-transition nodes (length(f_children) == 1 && f_deme[u] != f_deme[parent])
        jumps_12 = 0
        jumps_21 = 0
        jtimes_12 = Float64[]
        jtimes_21 = Float64[]

        kept_internals = Int[]
        function collect_induced(u_id::Int, parent_deme::Int)
            d_u = f_deme[u_id]
            if parent_deme > 0 && d_u != parent_deme
                if parent_deme == 1 && d_u == 2
                    jumps_12 += 1
                    push!(jtimes_12, f_time[u_id])
                elseif parent_deme == 2 && d_u == 1
                    jumps_21 += 1
                    push!(jtimes_21, f_time[u_id])
                end
            end
            if !is_sampled_leaf[u_id]
                ch = f_children[u_id]
                if length(ch) >= 2 || (parent_deme > 0 && d_u != parent_deme)
                    push!(kept_internals, u_id)
                    for c_id in ch
                        collect_induced(c_id, d_u)
                    end
                elseif length(ch) == 1
                    collect_induced(ch[1], parent_deme > 0 ? parent_deme : d_u)
                end
            end
        end
        collect_induced(orig_root, 0)

        n_int = length(kept_internals)
        n_tot = n_tips + n_int
        old_to_new = Dict{Int,Int}()
        new_node_state = zeros(Int, n_tot)
        new_node_time = zeros(Float64, n_tot)

        for (i, lf) in enumerate(sampled_leaf_nodes)
            old_to_new[lf] = i
            new_node_state[i] = f_deme[lf]
            new_node_time[i] = f_time[lf]
        end
        for (j, u_id) in enumerate(kept_internals)
            new_u = n_tips + j
            old_to_new[u_id] = new_u
            new_node_state[new_u] = f_deme[u_id]
            new_node_time[new_u] = f_time[u_id]
        end

        edges_p = Int[]
        edges_c = Int[]
        edge_len = Float64[]
        for u_id in kept_internals
            p_new = old_to_new[u_id]
            for c0 in f_children[u_id]
                curr = c0
                while !haskey(old_to_new, curr) && length(f_children[curr]) == 1
                    curr = f_children[curr][1]
                end
                c_new = old_to_new[curr]
                bl = max(1e-6, new_node_time[c_new] - new_node_time[p_new])
                push!(edges_p, p_new)
                push!(edges_c, c_new)
                push!(edge_len, bl)
            end
        end

        full_tree = PhyloTree(
            hcat(edges_p, edges_c), edge_len, tip_label, tip_state,
            new_node_state, new_node_time, n_int,
            jumps_12, jumps_21, sort!(jtimes_12), sort!(jtimes_21)
        )
        binary_tree = collapse_singletons(full_tree)
        return (
            tree = remove_singletons ? binary_tree : full_tree,
            full_tree = full_tree,
            binary_tree = binary_tree,
            final_state = copy(u),
            attempt = attempt
        )
    end
    error("simulate_forward_master_tree did not reach min_final_infected in $max_attempts attempts")
end

# ==============================================================================
# ModelingToolkit / SciML-style TreeProblem & Steppers
# ==============================================================================

"""
    TreeProblem(rs::ReactionSystem, infected_species::AbstractVector, u0, tspan, p, samples)

SciML / `ModelingToolkit.jl`-style problem specification for simulating phylogenetic trees
from a `ReactionSystem` using either:
- `CoalescentStepper()` (deterministic ODE + Volz Markovian structured coalescent, like `Coalescent.jl`)
- `ForwardSSAStepper()` (forward-in-time stochastic simulation tracking individual live lineages, like MASTER)
- `ReMASTERStepper()` (forward stochastic reaction log + backward hypergeometric event tracing, like ReMASTER)
"""
struct TreeProblem
    rs::ReactionSystem
    infected::Vector{Any}
    u0::Any
    tspan::Tuple{Float64,Float64}
    p::Any
    samples::SampleConfiguration
end

function TreeProblem(
    rs::ReactionSystem,
    infected::AbstractVector,
    u0,
    tspan::Tuple{<:Real,<:Real},
    p,
    samples
)
    scfg = samples isa SampleConfiguration ? samples : SampleConfiguration(samples, tspan[2])
    return TreeProblem(rs, collect(Any, infected), u0, (Float64(tspan[1]), Float64(tspan[2])), p, scfg)
end

"""`CoalescentStepper()`: solve a `TreeProblem` by deterministic demography and backward structured coalescence."""
struct CoalescentStepper end
"""`ForwardSSAStepper(; remove_singletons=true)`: solve a `TreeProblem` by tracking individual lineages in a forward Gillespie simulation."""
struct ForwardSSAStepper
    remove_singletons::Bool
end
ForwardSSAStepper(; remove_singletons::Bool=true) = ForwardSSAStepper(remove_singletons)

"""`ReMASTERStepper(; remove_singletons=true)`: solve a `TreeProblem` by forward reaction logging and backward event tracing."""
struct ReMASTERStepper
    remove_singletons::Bool
end
ReMASTERStepper(; remove_singletons::Bool=true) = ReMASTERStepper(remove_singletons)

function SciMLBase.solve(prob::TreeProblem, ::CoalescentStepper; rng::AbstractRNG=Random.default_rng(), solver=Tsit5())
    mod = ModelFGY(prob.rs, prob.infected, prob.u0, prob.tspan, prob.p)
    st = SimTree(mod, prob.samples; computedescendants=true, solver=solver, rng=rng)
    return simtree_to_phylotree(st; T_sample=prob.tspan[2], demes=mod.demes)
end

function SciMLBase.solve(prob::TreeProblem, stepper::ForwardSSAStepper; rng::AbstractRNG=Random.default_rng())
    res = simulate_forward_master_tree(
        prob.rs, prob.infected, prob.u0, prob.tspan, prob.p, prob.samples;
        remove_singletons=stepper.remove_singletons, rng=rng
    )
    return res.tree
end

function SciMLBase.solve(prob::TreeProblem, stepper::ReMASTERStepper; rng::AbstractRNG=Random.default_rng())
    res = simulate_remaster_tree(prob.rs, prob.infected, prob.u0, prob.tspan,
        prob.p, prob.samples; remove_singletons=stepper.remove_singletons, rng=rng)
    return res.tree
end

"""
    SimTree(rs::ReactionSystem, infected::AbstractVector, u0, tspan::Tuple{<:Real,<:Real}, p, samples; method::Symbol=:coalescent, rng=Random.default_rng())

Simulate a tree directly from a `ModelingToolkit.ReactionSystem` using ModelingToolkit-style
`u0`, `tspan`, `p`, and `samples` (`[Y1 => 100, Y2 => 100]`):
- `method = :coalescent` (default): deterministic ODE + Volz structured coalescent (like `Coalescent.jl`)
- `method = :forward`: forward-in-time stochastic simulation tracking individual live lineages (like MASTER)
- `method = :remaster`: forward stochastic reaction log + backward hypergeometric event tracing (like ReMASTER)
"""
function SimTree(
    rs::ReactionSystem,
    infected::AbstractVector,
    u0,
    tspan::Tuple{<:Real,<:Real},
    p,
    samples;
    method::Symbol = :coalescent,
    remove_singletons::Bool = true,
    rng::AbstractRNG = Random.default_rng()
)
    tprob = TreeProblem(rs, infected, u0, tspan, p, samples)
    if method === :coalescent
        return solve(tprob, CoalescentStepper(); rng=rng)
    elseif method === :forward || method === :master
        return solve(tprob, ForwardSSAStepper(remove_singletons); rng=rng)
    elseif method === :remaster || method === :stochastic || method === :forward_backward
        return solve(tprob, ReMASTERStepper(remove_singletons); rng=rng)
    else
        throw(ArgumentError("Unknown tree simulation method: $method"))
    end
end

"""Alias for [`ReMASTERStepper`](@ref), the forward-backward SSA tree solver."""
const ForwardBackwardSSAStepper = ReMASTERStepper
"""Alias for [`simulate_forward_master_tree`](@ref); accepts the same arguments and keywords."""
const simulate_forward_lineage_tree = simulate_forward_master_tree
