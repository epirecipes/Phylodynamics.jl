const _SANKOFF_INF = typemax(Int)

struct _ParsimonyAutomaton
    costs::Matrix{Int}
    messages::Vector{Vector{Int}}
    leaf_classes::Vector{Int}
    parent::Matrix{Int}
    increments::Matrix{Int}
end

_sankoff_add(a::Int, b::Int) =
    a == _SANKOFF_INF || b == _SANKOFF_INF ? _SANKOFF_INF : Base.Checked.checked_add(a, b)

function _sankoff_merge(costs, a, b)
    combined = _sankoff_add.(a, b)
    increment = minimum(combined)
    increment == _SANKOFF_INF &&
        throw(ArgumentError("parsimony costs admit a subtree with no possible ancestral state"))
    message = [minimum(_sankoff_add(costs[i, j], combined[j]) for j in eachindex(a))
               for i in eachindex(a)]
    return [v == _SANKOFF_INF ? v : v - increment for v in message], increment
end

function _compile_sankoff(costs::AbstractMatrix; max_classes::Int)
    k = size(costs, 1)
    k > 0 && size(costs, 2) == k ||
        throw(ArgumentError("parsimony costs must be a nonempty square matrix"))
    C = Matrix{Int}(undef, k, k)
    for i in 1:k, j in 1:k
        x = costs[i, j]
        x isa Real && (x == Inf || (isfinite(x) && x >= 0 && isinteger(x))) ||
            throw(ArgumentError("parsimony costs must be nonnegative integers or Inf"))
        C[i, j] = x == Inf ? _SANKOFF_INF : Int(x)
        x != Inf && C[i, j] == _SANKOFF_INF &&
            throw(ArgumentError("finite parsimony costs must be smaller than typemax(Int)"))
    end
    all(iszero, diag(C)) || throw(ArgumentError("parsimony costs must have zero diagonal"))
    max_classes > 0 || throw(ArgumentError("parsimony class budget must be positive"))

    messages = Vector{Int}[]
    indices = Dict{Tuple{Vararg{Int}},Int}()
    function register!(message)
        key = Tuple(message)
        haskey(indices, key) && return indices[key]
        length(messages) < max_classes ||
            throw(ArgumentError("parsimony closure exceeds $max_classes classes; increase " *
                "max_parsimony_classes/max_states, group character_states, or omit this group"))
        push!(messages, message)
        indices[key] = length(messages)
        return length(messages)
    end
    leaves = [register!(collect(C[:, j])) for j in 1:k]
    a = 1
    while a <= length(messages)
        for b in 1:a
            message, _ = _sankoff_merge(C, messages[a], messages[b])
            register!(message)
        end
        a += 1
    end
    n = length(messages)
    parent = zeros(Int, n, n)
    increments = zeros(Int, n, n)
    for a in 1:n, b in 1:a
        message, increment = _sankoff_merge(C, messages[a], messages[b])
        parent[a, b] = parent[b, a] = indices[Tuple(message)]
        increments[a, b] = increments[b, a] = increment
    end
    return _ParsimonyAutomaton(C, messages, leaves, parent, increments)
end

struct _CoalescentMomentLayout
    m::Int
    k::Int
    groups::Dict{Symbol,UnitRange{Int}}
    programs::Dict{Symbol,_ParsimonyAutomaton}
    character_states::Vector{Int}
    cherry_pairs::Vector{Tuple{Int,Int}}
    nstates::Int
end

function _coalescent_moment_layout(m, states, character_states, costs;
    max_states::Int, max_parsimony_classes::Int)
    selected = states isa Symbol ? [states] : collect(states)
    allowed = (:lineages, :jumps, :coalescences, :branch_lengths, :cherries,
        :sackin, :fitch, :camin_sokal, :sankoff, :parsimony_lengths)
    all(s -> s in allowed, selected) ||
        throw(ArgumentError("unknown state group; choose from $allowed"))
    length(unique(selected)) == length(selected) ||
        throw(ArgumentError("state groups must not be repeated"))
    length(character_states) == m && all(x -> x isa Integer && x > 0, character_states) ||
        throw(ArgumentError("character_states must give a positive integer label for each deme"))
    labels = Int.(character_states)
    k = maximum(labels)
    sort(unique(labels)) == collect(1:k) ||
        throw(ArgumentError("character_states must use contiguous labels 1:k"))
    max_states >= m || throw(ArgumentError("max_states cannot hold the $m lineage states"))
    max_parsimony_classes > 0 || throw(ArgumentError("max_parsimony_classes must be positive"))
    groups = Dict{Symbol,UnitRange{Int}}()
    programs = Dict{Symbol,_ParsimonyAutomaton}()
    nstates = 0
    function allocate!(name, n)
        n <= max_states - nstates ||
            throw(ArgumentError("selected ODE groups exceed max_states=$max_states"))
        groups[name] = (nstates + 1):(nstates + n)
        nstates += n
    end
    allocate!(:lineages, m)
    pairs = Tuple{Int,Int}[]
    for group in (:jumps, :coalescences, :branch_lengths, :cherries, :sackin)
        group in selected || continue
        if group === :cherries
            allocate!(:external, m * k)
            allocate!(:cherries, k * (k + 1) ÷ 2)
            pairs = [(a, b) for a in 1:k for b in a:k]
        elseif group === :sackin
            allocate!(:descendants, m)
            allocate!(:sackin, 1)
        else
            allocate!(group, group === :branch_lengths ? m : m * m)
        end
    end
    requests = Pair{Symbol,Matrix}[]
    if :fitch in selected
        k <= max_parsimony_classes && k <= (max_states - nstates - 1) ÷ m ||
            throw(ArgumentError("Fitch singleton classes already exceed the state budget"))
        push!(requests, :fitch => ones(Int, k, k) - Matrix{Int}(I, k, k))
    end
    if :camin_sokal in selected
        k == 2 || throw(ArgumentError("camin_sokal requires two character states; " *
            "use sankoff with explicit directed costs for more states"))
        push!(requests, :camin_sokal_12 => [0.0 1.0; Inf 0.0])
        push!(requests, :camin_sokal_21 => [0.0 Inf; 1.0 0.0])
    end
    if :sankoff in selected
        costs isa AbstractMatrix && size(costs) == (k, k) ||
            throw(ArgumentError("sankoff requires a k-by-k costs matrix"))
        push!(requests, :sankoff => Matrix(costs))
    elseif costs !== nothing
        throw(ArgumentError("costs was supplied but :sankoff was not selected"))
    end
    track_lengths = :parsimony_lengths in selected
    track_lengths && isempty(requests) &&
        throw(ArgumentError("parsimony_lengths requires a selected parsimony group"))
    for (name, C) in requests
        budget = min(max_parsimony_classes, (max_states - nstates - 1) ÷ (m + track_lengths))
        program = _compile_sankoff(C; max_classes=budget)
        programs[name] = program
        allocate!(name, m * length(program.messages))
        allocate!(Symbol(name, :_score), 1)
        track_lengths && allocate!(Symbol(name, :_lengths), length(program.messages))
    end
    return _CoalescentMomentLayout(m, k, groups, programs, labels, pairs, nstates)
end

function _coalescent_moment_initial(layout, samples)
    u = zeros(Float64, layout.nstates)
    u[layout.groups[:lineages]] .= samples
    if haskey(layout.groups, :external)
        L = reshape(view(u, layout.groups[:external]), layout.m, layout.k)
        for i in 1:layout.m
            L[i, layout.character_states[i]] = samples[i]
        end
    end
    haskey(layout.groups, :descendants) &&
        (u[layout.groups[:descendants]] .= samples)
    for (name, program) in layout.programs
        H = reshape(view(u, layout.groups[name]), layout.m, length(program.messages))
        for i in 1:layout.m
            H[i, program.leaf_classes[layout.character_states[i]]] = samples[i]
        end
    end
    return u
end

function _coalescent_moment_rhs!(du, u, F, G, Y, layout;
    descendant_mean = (x, a) -> iszero(a) ? zero(x) : x / a)
    fill!(du, 0)
    m, k, groups = layout.m, layout.k, layout.groups
    A = view(u, groups[:lineages])
    dA = view(du, groups[:lineages])
    p = A ./ Y
    for i in 1:m, j in 1:m
        coal = F[i, j] * p[i] * p[j]
        transfer = (F[i, j] * (1 - p[i]) + G[i, j]) * p[j]
        dA[j] -= coal
        if i != j
            dA[i] += transfer
            dA[j] -= transfer
        end
        if haskey(groups, :jumps) && i != j
            du[groups[:jumps][i + (j - 1) * m]] = (F[i, j] + G[i, j]) * p[j]
        end
        if haskey(groups, :coalescences)
            du[groups[:coalescences][i + (j - 1) * m]] = coal
        end
    end
    if haskey(groups, :branch_lengths)
        du[groups[:branch_lengths]] .= A
    end
    if haskey(groups, :external)
        L = reshape(view(u, groups[:external]), m, k)
        dL = reshape(view(du, groups[:external]), m, k)
        for i in 1:m, j in 1:m
            for a in 1:k
                dL[i, a] -= F[i, j] * (L[i, a] / Y[i]) * p[j]
                dL[j, a] -= F[i, j] * p[i] * (L[j, a] / Y[j])
                if i != j
                    transfer = (F[i, j] * (1 - p[i]) + G[i, j]) * L[j, a] / Y[j]
                    dL[i, a] += transfer
                    dL[j, a] -= transfer
                end
            end
            for (index, (a, b)) in enumerate(layout.cherry_pairs)
                rate = L[i, a] * L[j, b]
                a != b && (rate += L[i, b] * L[j, a])
                du[groups[:cherries][index]] += F[i, j] * rate / (Y[i] * Y[j])
            end
        end
    end
    if haskey(groups, :descendants)
        X = view(u, groups[:descendants])
        dX = view(du, groups[:descendants])
        means = descendant_mean.(X, A)
        for i in 1:m, j in 1:m
            if i != j
                transfer = (F[i, j] + G[i, j]) * p[j] * means[j]
                dX[i] += transfer
                dX[j] -= transfer
            end
            du[only(groups[:sackin])] += F[i, j] * p[i] * p[j] * (means[i] + means[j])
        end
    end
    for (name, program) in layout.programs
        n = length(program.messages)
        H = reshape(view(u, groups[name]), m, n)
        dH = reshape(view(du, groups[name]), m, n)
        score = only(groups[Symbol(name, :_score)])
        lengths = Symbol(name, :_lengths)
        if haskey(groups, lengths)
            du[groups[lengths]] .= vec(sum(H; dims=1))
        end
        for i in 1:m, j in 1:m
            if i != j
                for b in 1:n
                    transfer = (F[i, j] * (1 - p[i]) + G[i, j]) * H[j, b] / Y[j]
                    dH[i, b] += transfer
                    dH[j, b] -= transfer
                end
            end
            for a in 1:n, b in 1:n
                rate = F[i, j] * (H[i, a] / Y[i]) * (H[j, b] / Y[j])
                dH[i, a] -= rate
                dH[j, b] -= rate
                dH[i, program.parent[a, b]] += rate
                du[score] += program.increments[a, b] * rate
            end
        end
    end
    return nothing
end

function _validate_moment_tfgy(tfgy)
    m, n = length(tfgy.demes), length(tfgy.times)
    m > 0 && length(unique(tfgy.demes)) == m ||
        throw(ArgumentError("TFGY must have nonempty, unique deme names"))
    n >= 2 && all(isfinite, tfgy.times) && all(diff(tfgy.times) .< 0) ||
        throw(ArgumentError("TFGY times must be finite and strictly decreasing, with at least two points"))
    all(v -> length(v) == n, (tfgy.Fs, tfgy.Gs, tfgy.Ys, tfgy.deaths)) ||
        throw(ArgumentError("TFGY arrays must have the same length as times"))
    for h in 1:n
        size(tfgy.Fs[h]) == (m, m) && size(tfgy.Gs[h]) == (m, m) &&
            length(tfgy.Ys[h]) == m && length(tfgy.deaths[h]) == m ||
            throw(ArgumentError("TFGY dimensions do not match its deme count"))
        all(x -> isfinite(x) && x >= 0, tfgy.Fs[h]) &&
            all(x -> isfinite(x) && x >= 0, tfgy.Gs[h]) &&
            all(x -> isfinite(x) && x >= 0, tfgy.deaths[h]) ||
            throw(ArgumentError("TFGY rates must be finite and nonnegative"))
        all(x -> isfinite(x) && x >= 0, tfgy.Ys[h]) ||
            throw(ArgumentError("TFGY population sizes must be finite and nonnegative"))
    end
    return nothing
end

function _moment_fgy_at(tfgy, t)
    h = clamp(searchsortedlast(tfgy.times, t; rev=true), 1, length(tfgy.times) - 1)
    w = (t - tfgy.times[h + 1]) / (tfgy.times[h] - tfgy.times[h + 1])
    return (w .* tfgy.Fs[h] .+ (1 - w) .* tfgy.Fs[h + 1],
        w .* tfgy.Gs[h] .+ (1 - w) .* tfgy.Gs[h + 1],
        w .* tfgy.Ys[h] .+ (1 - w) .* tfgy.Ys[h + 1])
end

function _moment_cherry_matrix(values, layout)
    C = zeros(eltype(values), layout.k, layout.k)
    for (index, (a, b)) in enumerate(layout.cherry_pairs)
        C[a, b] = C[b, a] = values[index]
    end
    return C
end

"""
    solve_coalescent_odes(tfgy::TFGY, samples; states=(:lineages,), kwargs...)
    solve_coalescent_odes(rs::ReactionSystem, infected, u0, tspan, p, samples; res=1000, kwargs...)

Integrate multi-deme, homochronous mean-field coalescent moment equations backwards
from `T_sample` (default `first(tfgy.times)`) to `t0` (default `last(tfgy.times)`),
stopping at total lineage count one if it is reached. `TFGY` rates are linearly
interpolated; its population sizes must be strictly positive throughout.

`states` selects ODE groups, not just output columns. Lineages are always included.
Optional groups are `:jumps`, `:coalescences`, `:branch_lengths`, `:cherries`,
`:sackin`, `:fitch`, `:camin_sokal` (both binary directional priors), and `:sankoff`.
`:parsimony_lengths` additionally integrates lengths by class for each selected
parsimony program (summed over demographic demes).
Cherries allocate external-lineage counts; Sackin allocates descendant masses.
`character_states` maps each deme to a contiguous character label `1:k` (default
`1:m`), allowing fewer character states than demes.

`:sankoff` requires a `costs` matrix of nonnegative integer directed edge costs
with zero diagonal; `Inf` forbids an edge. A free root is assumed. The reachable
normalized cost-vector closure is constructed exactly, or fails before integration
if incompatible or above budget. `max_states=10_000` bounds total ODE variables;
`max_parsimony_classes=256` bounds each parsimony automaton. Fitch can require
`2^k-1` classes per deme and parsimony updates are quadratic in the class count.
No parsimony classes are built unless requested.

Returns `(summary, traj, sol, layout, fwd)`. `layout.groups` indexes the selected
ODE blocks; `layout.programs` records cost messages, leaf classes and merge tables.
`summary.reached_mrca` is false if only a partial ancestral history was covered.
`traj` contains only the allocated states, with columns `group_index`. Matrices in
`summary` use parent/source rows and child/destination columns; cherry matrices
store each unordered off-diagonal count in both symmetric entries.

These are deterministic moment-closure equations, not exact stochastic expectations.
The `Phylodynamics.Zoonosis` submodule supplies the study-specific
`ZoonosisParams` specialization with its two-host result schema.
"""
function solve_coalescent_odes(tfgy::TFGY, samples::AbstractVector{<:Real};
    states = (:lineages,), character_states = collect(1:length(tfgy.demes)),
    costs = nothing, max_states::Int = 10_000, max_parsimony_classes::Int = 256,
    T_sample::Union{Nothing,Real} = nothing, t0::Union{Nothing,Real} = nothing,
    dt::Real = 0.01, solver = Tsit5(), reltol::Real = 1e-9, abstol::Real = 1e-11)
    _validate_moment_tfgy(tfgy)
    T_sample === nothing && (T_sample = first(tfgy.times))
    t0 === nothing && (t0 = last(tfgy.times))
    m = length(tfgy.demes)
    length(samples) == m && all(x -> isfinite(x) && x >= 0, samples) && sum(samples) > 1 ||
        throw(ArgumentError("samples must be nonnegative deme counts with total greater than one"))
    isfinite(T_sample) && isfinite(t0) && last(tfgy.times) <= t0 < T_sample <= first(tfgy.times) ||
        throw(ArgumentError("sampling/origin times must lie inside the TFGY interval with t0 < T_sample"))
    isfinite(dt) && dt > 0 && isfinite(reltol) && reltol > 0 && isfinite(abstol) && abstol > 0 ||
        throw(ArgumentError("dt and solver tolerances must be finite and positive"))
    _, _, Y_sample = _moment_fgy_at(tfgy, T_sample)
    _, _, Y_origin = _moment_fgy_at(tfgy, t0)
    all(>(0), Y_sample) && all(>(0), Y_origin) &&
        all(h -> !(t0 <= tfgy.times[h] <= T_sample) || all(>(0), tfgy.Ys[h]),
            eachindex(tfgy.times)) ||
        throw(ArgumentError("analytical moment ODEs require strictly positive Y throughout the selected interval"))
    all(samples .<= Y_sample) || throw(ArgumentError("sample counts exceed population sizes at T_sample"))
    layout = _coalescent_moment_layout(m, states, character_states, costs;
        max_states, max_parsimony_classes)
    u0 = _coalescent_moment_initial(layout, samples)
    function rhs!(du, u, p, s)
        F, G, Y = _moment_fgy_at(tfgy, clamp(T_sample - s, t0, T_sample))
        _coalescent_moment_rhs!(du, u, F, G, Y, layout)
    end
    duration = Float64(T_sample - t0)
    grid = collect(0.0:Float64(dt):duration)
    last(grid) < duration && push!(grid, duration)
    root = ContinuousCallback((u, t, integrator) -> sum(view(u, layout.groups[:lineages])) - 1,
        nothing, integrator -> terminate!(integrator); save_positions=(false, true))
    problem = ODEProblem(rhs!, u0, (0.0, duration))
    sol = solve(problem, solver; saveat=grid, callback=root, reltol, abstol)
    SciMLBase.successful_retcode(sol) || error("coalescent moment integration failed: $(sol.retcode)")
    for (s, u) in zip(sol.t, sol.u)
        _, _, Y = _moment_fgy_at(tfgy, clamp(T_sample - s, t0, T_sample))
        A = view(u, layout.groups[:lineages])
        tolerance = 100 .* (abstol .+ reltol .* max.(Y, 1))
        all((-tolerance .<= A) .& (A .<= Y .+ tolerance)) ||
            throw(DomainError(T_sample - s,
                "moment approximation left 0 <= A <= Y; restrict the interval or revise the sampling/demography"))
    end
    final = last(sol.u)
    summary = Dict{Symbol,Any}(
        :T_sample => Float64(T_sample), :t_end => Float64(T_sample - last(sol.t)),
        :reached_mrca => sol.retcode == SciMLBase.ReturnCode.Terminated,
        :demes => copy(tfgy.demes), :samples => Float64.(samples))
    traj = DataFrame(s=sol.t, t_fwd=T_sample .- sol.t)
    for (name, range) in sort!(collect(layout.groups); by=p -> first(last(p)))
        values = final[range]
        if name in (:jumps, :coalescences)
            summary[name] = reshape(values, m, m)
        elseif name === :external
            summary[name] = reshape(values, m, layout.k)
        elseif haskey(layout.programs, name)
            summary[name] = reshape(values, m, length(layout.programs[name].messages))
        elseif name === :cherries
            summary[name] = _moment_cherry_matrix(values, layout)
        elseif name === :sackin || endswith(string(name), "_score")
            summary[name] = only(values)
        else
            summary[name] = values
        end
        for (index, state) in enumerate(range)
            traj[!, Symbol(name, :_, index)] = [u[state] for u in sol.u]
        end
    end
    return (summary=(; summary...), traj=traj, sol=sol, layout=layout, fwd=tfgy)
end

function solve_coalescent_odes(rs::ReactionSystem, infected::AbstractVector,
    u0, tspan::Tuple{<:Real,<:Real}, p, samples::AbstractVector; res::Int=1000, kwargs...)
    res >= 2 || throw(ArgumentError("res must be at least two"))
    all(isfinite, tspan) && tspan[1] < tspan[2] ||
        throw(ArgumentError("tspan must be a finite increasing interval"))
    counts = if eltype(samples) <: Pair
        names = [_species_name(s) for s in infected]
        mapped = Dict(_species_name(first(pair)) => last(pair) for pair in samples)
        length(mapped) == length(samples) && all(name -> name in names, keys(mapped)) ||
            throw(ArgumentError("sample map contains repeated or unknown infected species"))
        [get(mapped, name, 0.0) for name in names]
    else
        samples
    end
    tfgy = reaction_system_to_fgy(rs, infected, u0, p, Float64.(tspan); res)
    first(tfgy.times) == tspan[2] && last(tfgy.times) == tspan[1] ||
        error("forward demography did not cover the requested time interval")
    return solve_coalescent_odes(tfgy, counts; kwargs...)
end
