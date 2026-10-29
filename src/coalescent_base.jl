# ==============================================================================
# coalescent_base.jl
# Full implementation and modernization of Erik Volz's Coalescent.jl
# (https://github.com/emvolz/Coalescent.jl):
# - Reaction types: RXN_BIRTH, RXN_MIG, RXN_DEATH, RXN_DYNVAR
# - Coalescent events: SAMPLE, COALESCENT, MIGRATION, RECOMBINATION
# - Algorithms: ALGO_MARKOV, ALGO_STATIONARY
# - SampleConfiguration, CoalescentReaction, ModelFGY, solveodes
# - Event, SimTree, tonewick, distancematrix
# - Bidirectional conversion between SimTree and PhyloTree
# ==============================================================================

"""Reaction code for transmission/birth events in a `ModelFGY`."""
const RXN_BIRTH = 0
"""Reaction code for physical migration between demes."""
const RXN_MIG = 1
"""Reaction code for infected-lineage removal."""
const RXN_DEATH = 2
"""Reaction code for changes to auxiliary demographic variables."""
const RXN_DYNVAR = 3

"""`Event` code for a sampled lineage."""
const SAMPLE = 0
"""`Event` code for a coalescence."""
const COALESCENT = 1
"""`Event` code for a lineage changing demes."""
const MIGRATION = 2
"""`Event` code reserved for recombination."""
const RECOMBINATION = 3

"""Markovian lineage-state approximation for `SimTree`."""
const ALGO_MARKOV = "markov"
"""Stationary lineage-state approximation for `SimTree`."""
const ALGO_STATIONARY = "stationary"

"""
    SampleConfiguration

Represents a configuration for sampling in coalescent simulations (from `Coalescent.jl`).

# Fields
- `sconf::Vector{Tuple{Union{Nothing,String}, Float64}}`: Vector of `(deme, sampling_time)` tuples.
"""
struct SampleConfiguration
    sconf::Vector{Tuple{Union{Nothing,String}, Float64}}
end

function SampleConfiguration(demes::AbstractVector{<:AbstractString}, times::AbstractVector{<:Real})
    length(demes) == length(times) || throw(ArgumentError("demes and times must have the same length"))
    return SampleConfiguration([(String(d), Float64(t)) for (d, t) in zip(demes, times)])
end

function SampleConfiguration(deme_counts::Dict{String,Int}, t_sample::Real)
    sconf = Tuple{Union{Nothing,String}, Float64}[]
    for (d, n) in sort!(collect(deme_counts); by=first)
        for _ in 1:n
            push!(sconf, (d, Float64(t_sample)))
        end
    end
    return SampleConfiguration(sconf)
end

_species_name(x) = x isa AbstractString ? String(x) :
                   x isa Symbol ? String(x) :
                   String(Symbolics.tosymbol(x; escape=false))

"""
    SampleConfiguration(samples::AbstractVector{<:Pair}, t_default::Real=0.0)

Construct a `SampleConfiguration` using ModelingToolkit-style `Pair` syntax:
- `[Y1 => 100, Y2 => 100]` with `t_default = 12.85` (homochronous sampling at `t_default`)
- `[Y1 => (100, 12.85), Y2 => (100, 12.85)]` (count and sampling time tuple)
- `[Y1 => fill(12.85, 100), Y2 => fill(12.85, 100)]` (explicit vector of sample times per deme)
"""
function SampleConfiguration(samples::AbstractVector{<:Pair}, t_default::Real = 0.0)
    sconf = Tuple{Union{Nothing,String}, Float64}[]
    for (sp, val) in samples
        dname = _species_name(sp)
        if val isa Integer
            for _ in 1:Int(val)
                push!(sconf, (dname, Float64(t_default)))
            end
        elseif val isa Tuple{<:Integer, <:Real}
            cnt, t_val = val
            for _ in 1:Int(cnt)
                push!(sconf, (dname, Float64(t_val)))
            end
        elseif val isa AbstractVector{<:Real}
            for t_val in val
                push!(sconf, (dname, Float64(t_val)))
            end
        else
            throw(ArgumentError("Unsupported sample specification for $dname: $val"))
        end
    end
    return SampleConfiguration(sconf)
end

"""
    CoalescentReaction

Represents a birth, migration, death, or non-deme dynamic variable reaction in `ModelFGY`
(named `CoalescentReaction` to avoid name collision with `ModelingToolkit.Reaction`).
"""
mutable struct CoalescentReaction
    source::String
    recipient::Union{Nothing,String}
    type::Int
    expr::Union{Symbol,Expr,Number}
    rate_fn::Union{Nothing,Function}

    function CoalescentReaction(
        s::AbstractString,
        r::Union{Nothing,AbstractString},
        t::Int,
        e::Union{Symbol,Expr,Number};
        rate_fn::Union{Nothing,Function}=nothing
    )
        @assert t in (RXN_BIRTH, RXN_MIG, RXN_DEATH, RXN_DYNVAR)
        return new(String(s), r === nothing ? nothing : String(r), t, e, rate_fn)
    end
end

function CoalescentReaction(s::AbstractString, t::Int, e::Union{Symbol,Expr,Number}; rate_fn=nothing)
    @assert t in (RXN_DEATH, RXN_DYNVAR)
    return CoalescentReaction(s, nothing, t, e; rate_fn=rate_fn)
end

"""
    ModelFGY

Structured Forward-in-time Genealogy (FGY) model for coalescent simulations and inference
(compatible with `Coalescent.jl`).
"""
mutable struct ModelFGY
    modelname::String
    birthrxn::Vector{CoalescentReaction}
    migrationrxn::Vector{CoalescentReaction}
    deathrxn::Vector{CoalescentReaction}
    nondemerxn::Vector{CoalescentReaction}
    demes::Vector{String}
    nondemes::Vector{String}
    numberdemes::Int
    numbernondemes::Int
    initial::Dict{String,Float64}
    t0::Float64
    tfin::Float64
    parameters::Dict{String,Float64}
    helperexprs::Union{Nothing,Vector{Expr}}
end

function ModelFGY(;
    modelname::AbstractString = "model",
    birthrxn::Vector{CoalescentReaction} = CoalescentReaction[],
    migrationrxn::Vector{CoalescentReaction} = CoalescentReaction[],
    deathrxn::Vector{CoalescentReaction} = CoalescentReaction[],
    nondemerxn::Vector{CoalescentReaction} = CoalescentReaction[],
    demes::Vector{String} = String[],
    nondemes::Vector{String} = String[],
    initial::Dict{String,Float64} = Dict{String,Float64}(),
    t0::Real = 0.0,
    tfin::Real = 12.85,
    parameters::Dict{String,Float64} = Dict{String,Float64}(),
    helperexprs::Union{Nothing,Vector{Expr}} = nothing
)
    if isempty(demes)
        dset = String[]
        for r in birthrxn
            !(r.source in dset) && push!(dset, r.source)
            r.recipient !== nothing && !(r.recipient in dset) && push!(dset, r.recipient)
        end
        for r in migrationrxn
            !(r.source in dset) && push!(dset, r.source)
            r.recipient !== nothing && !(r.recipient in dset) && push!(dset, r.recipient)
        end
        for r in deathrxn
            !(r.source in dset) && push!(dset, r.source)
        end
        demes = dset
    end
    if isempty(nondemes) && !isempty(initial)
        nondemes = [k for k in keys(initial) if !(k in demes)]
    end
    return ModelFGY(
        String(modelname), birthrxn, migrationrxn, deathrxn, nondemerxn,
        demes, nondemes, length(demes), length(nondemes),
        initial, Float64(t0), Float64(tfin), parameters, helperexprs
    )
end

function Base.show(io::IO, x::ModelFGY)
    print(io, "ModelFGY(\"$(x.modelname)\", demes=$(x.demes), nondemes=$(x.nondemes), t=($(x.t0), $(x.tfin)))")
end

"""
    ModelFGY(params::ZoonosisParams; modelname::String="zoonosis") -> ModelFGY

Construct a `Coalescent.jl` `ModelFGY` for the two-host zoonosis model (Models 1, 2, or 3),
where `demes = ["Y1", "Y2"]` and `nondemes = ["X1", "X2"]`.
"""
function ModelFGY(p::ZoonosisParams; modelname::String="zoonosis", t0::Real=0.0, tfin::Real=p.t_peak)
    brxns = [
        CoalescentReaction("Y1", "Y1", RXN_BIRTH, :(beta * X1 * Y1);
            rate_fn = (u, par, t) -> par[:beta] * max(u[3], 0.0) * max(u[1], 0.0)),
        CoalescentReaction("Y1", "Y2", RXN_BIRTH, :(epsilon * beta * X2 * Y1);
            rate_fn = (u, par, t) -> par[:epsilon] * par[:beta] * max(u[4], 0.0) * max(u[1], 0.0)),
        CoalescentReaction("Y2", "Y2", RXN_BIRTH, :(phi * beta * X2 * Y2);
            rate_fn = (u, par, t) -> par[:phi] * par[:beta] * max(u[4], 0.0) * max(u[2], 0.0)),
        CoalescentReaction("Y2", "Y1", RXN_BIRTH, :(psi * beta * X1 * Y2);
            rate_fn = (u, par, t) -> par[:psi] * par[:beta] * max(u[3], 0.0) * max(u[2], 0.0))
    ]
    migrxns = CoalescentReaction[]
    deathrxns = [
        CoalescentReaction("Y1", RXN_DEATH, :((gamma1 + mu1 + m1) * Y1);
            rate_fn = (u, par, t) -> (par[:gamma1] + par[:mu1] + par[:m1]) * max(u[1], 0.0)),
        CoalescentReaction("Y2", RXN_DEATH, :((gamma2 + mu2 + m2) * Y2);
            rate_fn = (u, par, t) -> (par[:gamma2] + par[:mu2] + par[:m2]) * max(u[2], 0.0))
    ]
    nondemerxns = [
        CoalescentReaction("X1", RXN_DYNVAR,
            :(nu1 - beta * X1 * Y1 - psi * beta * X1 * Y2 - mu1 * X1);
            rate_fn = (u, par, t) -> par[:nu1] - par[:beta] * u[3] * u[1] - par[:psi] * par[:beta] * u[3] * u[2] - par[:mu1] * u[3]),
        CoalescentReaction("X2", RXN_DYNVAR,
            :(nu2 - epsilon * beta * X2 * Y1 - phi * beta * X2 * Y2 - mu2 * X2);
            rate_fn = (u, par, t) -> par[:nu2] - par[:epsilon] * par[:beta] * u[4] * u[1] - par[:phi] * par[:beta] * u[4] * u[2] - par[:mu2] * u[4])
    ]
    demes = ["Y1", "Y2"]
    nondemes = ["X1", "X2"]
    initial = Dict(
        "Y1" => p.Y1_0,
        "Y2" => p.Y2_0,
        "X1" => p.X1_0,
        "X2" => p.X2_0
    )
    parameters = Dict(
        "nu1" => p.nu1, "mu1" => p.mu1, "gamma1" => p.gamma1, "m1" => p.m1,
        "beta" => p.beta, "nu2" => p.nu2, "mu2" => p.mu2, "gamma2" => p.gamma2,
        "m2" => p.m2, "epsilon" => p.epsilon, "phi" => p.phi, "psi" => p.psi
    )
    return ModelFGY(;
        modelname=modelname, birthrxn=brxns, migrationrxn=migrxns,
        deathrxn=deathrxns, nondemerxn=nondemerxns, demes=demes, nondemes=nondemes,
        initial=initial, t0=Float64(t0), tfin=Float64(tfin), parameters=parameters
    )
end

# Evaluate a CoalescentReaction rate given state vector u (ordered [demes..., nondemes...]), parameter Dict, and time t
function _eval_rxn_rate(
    rxn::CoalescentReaction,
    u::AbstractVector{<:Real},
    var_index::Dict{String,Int},
    sym_params::NamedTuple,
    str_params::Dict{String,Float64},
    t::Real
)
    if rxn.rate_fn !== nothing
        return Float64(rxn.rate_fn(u, sym_params, t))
    end
    # Fallback expression evaluator for user-supplied Expr
    return _eval_expr(rxn.expr, u, var_index, str_params, t)
end

function _eval_expr(ex, u::AbstractVector{<:Real}, var_index::Dict{String,Int}, params::Dict{String,Float64}, t::Real)
    if ex isa Real
        return Float64(ex)
    elseif ex isa Symbol
        s = String(ex)
        s == "t" && return Float64(t)
        haskey(var_index, s) && return Float64(u[var_index[s]])
        haskey(params, s) && return Float64(params[s])
        return Float64(getfield(Base, ex))
    elseif ex isa Expr && ex.head === :call
        op = ex.args[1]
        args = [_eval_expr(a, u, var_index, params, t) for a in ex.args[2:end]]
        f = op isa Symbol ? getfield(Base, op) : op
        return Float64(f(args...))
    else
        error("Unsupported expression in CoalescentReaction: $ex")
    end
end

"""
    solveodes(model::ModelFGY; solver=Tsit5(), res::Union{Missing,Int}=missing, reltol=1e-9, abstol=1e-11)

Solve the forward demographic ODEs defined by `model::ModelFGY` (compatible with `Coalescent.jl`).
State vector order is `vcat(model.demes, model.nondemes)`.
"""
function solveodes(
    model::ModelFGY;
    solver = Tsit5(),
    res::Union{Missing,Int} = missing,
    reltol::Float64 = 1e-9,
    abstol::Float64 = 1e-11
)
    all_vars = vcat(model.demes, model.nondemes)
    var_index = Dict(v => i for (i, v) in enumerate(all_vars))
    deme_index = Dict(d => i for (i, d) in enumerate(model.demes))
    nondeme_index = Dict(nd => model.numberdemes + i for (i, nd) in enumerate(model.nondemes))
    sym_params = NamedTuple{Tuple(Symbol.(keys(model.parameters)))}(Tuple(values(model.parameters)))
    str_params = model.parameters

    function fgy_odes!(du, u, p, t)
        fill!(du, 0.0)
        for rxn in model.birthrxn
            r = _eval_rxn_rate(rxn, u, var_index, sym_params, str_params, t)
            du[deme_index[rxn.recipient]] += r
        end
        for rxn in model.migrationrxn
            r = _eval_rxn_rate(rxn, u, var_index, sym_params, str_params, t)
            du[deme_index[rxn.recipient]] += r
            du[deme_index[rxn.source]] -= r
        end
        for rxn in model.deathrxn
            r = _eval_rxn_rate(rxn, u, var_index, sym_params, str_params, t)
            du[deme_index[rxn.source]] -= r
        end
        for rxn in model.nondemerxn
            r = _eval_rxn_rate(rxn, u, var_index, sym_params, str_params, t)
            du[nondeme_index[rxn.source]] += r
        end
        return nothing
    end

    u0 = [Float64(model.initial[k]) for k in all_vars]
    tspan = (model.t0, model.tfin)
    prob = ODEProblem(fgy_odes!, u0, tspan)
    if ismissing(res)
        return solve(prob, solver; reltol=reltol, abstol=abstol)
    else
        saveat = collect(range(model.t0, model.tfin; length=res))
        return solve(prob, solver; saveat=saveat, reltol=reltol, abstol=abstol)
    end
end

"""
    Event

Represents an event in a coalescent process (`SAMPLE`, `COALESCENT`, `MIGRATION`, `RECOMBINATION`).
"""
mutable struct Event
    type::Int
    height::Float64
    source::Union{Nothing,String}
    sink::Union{Nothing,String}
end

Event(t::Int, b::Real) = Event(t, Float64(b), nothing, nothing)

"""
    SimTree

Simulated coalescent tree data structure from `Coalescent.jl`.
"""
Base.@kwdef mutable struct SimTree
    parent::Vector{Union{Nothing,Int}}
    child::Vector{Union{Nothing,Int}}
    n::Int
    nNode::Int
    edgelength::Vector{Union{Nothing,Float64}}
    heights::Vector{Union{Nothing,Float64}}
    tiplabs::Vector{Union{Nothing,String}}
    shs::Vector{Union{Nothing,Float64}}
    descendants::Union{Nothing,BitMatrix} = nothing
    daughters::Union{Nothing,Vector{Tuple{Int,Int,Int,Int}}} = nothing
    demes::Union{Nothing,Vector{String}} = nothing
    true_jumps_12::Int = 0
    true_jumps_21::Int = 0
    true_jump_times_12::Vector{Float64} = Float64[]
    true_jump_times_21::Vector{Float64} = Float64[]
end

function Base.show(io::IO, x::SimTree)
    print(io, "SimTree(n=$(x.n) tips, nNode=$(x.nNode) internal nodes, J12=$(x.true_jumps_12), J21=$(x.true_jumps_21))")
end

"""
    SimTree(events::Vector{Event}, model::ModelFGY; computedescendants::Bool=false, rng=Random.default_rng()) -> SimTree

Construct a `SimTree` from an array of genealogical `Event`s and a `ModelFGY`
(adapted from `Coalescent.jl/src/s0.jl`, tracking exact directional deme jumps).
"""
function SimTree(
    events_in::Vector{Event},
    model::ModelFGY;
    computedescendants::Bool = false,
    rng::AbstractRNG = Random.default_rng()
)::SimTree
    events = deepcopy(events_in)
    ix = sortperm([e.height for e in events])
    events = events[ix]
    n = count(e -> e.type == SAMPLE, events)
    nNode = n - 1
    nCo = count(e -> e.type == COALESCENT, events)
    mrcaheight = isempty(events) ? 0.0 : events[end].height
    epsilonheight = 1e-5
    cotoadd = nNode - nCo
    for i in 1:cotoadd
        push!(events, Event(COALESCENT, mrcaheight + epsilonheight * i))
    end

    tiplabs = Union{Nothing,String}[]
    shs = Union{Nothing,Float64}[]
    sample_idx = 0
    for e in events
        if e.type == SAMPLE
            sample_idx += 1
            suffix = (e.source == "Y2" || e.source == "Y") ? "Y" : "I"
            push!(tiplabs, "t.$(e.source).$(sample_idx)_$(suffix)")
            push!(shs, e.height)
        end
    end

    nedges = n + nNode - 1
    parent = Vector{Union{Nothing,Int}}(nothing, nedges)
    child = Vector{Union{Nothing,Int}}(nothing, nedges)
    daughters = Vector{Tuple{Int,Int,Int,Int}}(undef, nNode)
    edgelength = Vector{Union{Nothing,Float64}}(nothing, nedges)
    heights = Vector{Union{Nothing,Float64}}(nothing, n + nNode)
    demes = Vector{String}(undef, n + nNode)

    extant = Dict(d => Int[] for d in model.demes)
    for v in values(extant)
        sizehint!(v, n + nNode)
    end

    a = n + nNode
    su = 1
    ie = 1
    desc = computedescendants ? falses(nNode, n) : nothing

    jumps_12 = 0
    jumps_21 = 0
    jtimes_12 = Float64[]
    jtimes_21 = Float64[]
    mst = model.tfin

    for e in events
        if e.type == SAMPLE
            heights[su] = e.height
            push!(extant[e.source], su)
            demes[su] = e.source
            su += 1
        elseif e.type == COALESCENT
            heights[a] = e.height
            if isnothing(e.source)
                e.source = model.demes[1]
                e.sink = !isempty(extant[model.demes[end]]) && length(extant[model.demes[1]]) < 2 ? model.demes[end] : model.demes[1]
            end
            e0 = extant[e.source]
            e1 = extant[e.sink]
            if isempty(e0) || isempty(e1) || (e.source == e.sink && length(e0) < 2)
                # Fallback if forced coalescence at boundary
                all_ext = reduce(vcat, values(extant))
                length(all_ext) < 2 && continue
                u, v = sample(rng, all_ext, 2; replace=false)
                for d in model.demes
                    filter!(x -> x != u && x != v, extant[d])
                end
                e.source = model.demes[1]
                e.sink = model.demes[1]
            else
                iu = rand(rng, 1:length(e0))
                if e.source == e.sink
                    choices = [idx for idx in 1:length(e1) if idx != iu]
                    iv = rand(rng, choices)
                else
                    iv = rand(rng, 1:length(e1))
                end
                u = e0[iu]
                v = e1[iv]
                if e.source == e.sink
                    deleteat!(e0, sort!([iu, iv]))
                else
                    deleteat!(e0, iu)
                    deleteat!(e1, iv)
                end
            end

            # Record cross-deme transmission jump at coalescence if source != sink
            if e.source != e.sink
                t_fwd = max(0.0, mst - e.height)
                if e.source == "Y1" && e.sink == "Y2"
                    jumps_12 += 1
                    push!(jtimes_12, t_fwd)
                elseif e.source == "Y2" && e.sink == "Y1"
                    jumps_21 += 1
                    push!(jtimes_21, t_fwd)
                end
            end

            push!(extant[e.source], a)
            demes[a] = e.source

            parent[ie] = a
            child[ie] = u
            edgelength[ie] = heights[a] - heights[u]
            ie += 1
            parent[ie] = a
            child[ie] = v
            edgelength[ie] = heights[a] - heights[v]
            ie += 1
            daughters[a - n] = (u, v, ie - 2, ie - 1)

            if computedescendants
                if u <= n
                    desc[a - n, u] = true
                else
                    desc[a - n, :] .|= desc[u - n, :]
                end
                if v <= n
                    desc[a - n, v] = true
                else
                    desc[a - n, :] .|= desc[v - n, :]
                end
            end
            a -= 1
        elseif e.type == MIGRATION
            e1 = extant[e.sink]
            if !isempty(e1)
                iv = rand(rng, 1:length(e1))
                v = e1[iv]
                deleteat!(e1, iv)
                push!(extant[e.source], v)
                t_fwd = max(0.0, mst - e.height)
                if e.source == "Y1" && e.sink == "Y2"
                    jumps_12 += 1
                    push!(jtimes_12, t_fwd)
                elseif e.source == "Y2" && e.sink == "Y1"
                    jumps_21 += 1
                    push!(jtimes_21, t_fwd)
                end
            end
        end
    end

    return SimTree(;
        parent=parent, child=child, n=n, nNode=nNode,
        edgelength=edgelength, heights=heights, tiplabs=tiplabs, shs=shs,
        descendants=desc, daughters=daughters, demes=demes,
        true_jumps_12=jumps_12, true_jumps_21=jumps_21,
        true_jump_times_12=sort!(jtimes_12), true_jump_times_21=sort!(jtimes_21)
    )
end

"""
    SimTree(model::ModelFGY, sample::SampleConfiguration; computedescendants::Bool=false, solver=Tsit5(), ytol::Float64=1e-6, rng=Random.default_rng()) -> SimTree

Simulate a structured coalescent tree from `model::ModelFGY` and `sample::SampleConfiguration`
using the Markovian cumulative-hazard ODE + `ContinuousCallback` algorithm of `Coalescent.jl`.
"""
function SimTree(
    model::ModelFGY,
    sample_cfg::SampleConfiguration;
    computedescendants::Bool = false,
    solver = Tsit5(),
    ytol::Float64 = 1e-6,
    rng::AbstractRNG = Random.default_rng()
)::SimTree
    sampletimes = [Float64(x[2]) for x in sample_cfg.sconf]
    samplestates = [String(x[1]) for x in sample_cfg.sconf]
    n = length(sampletimes)
    @assert n > 1

    mst = maximum(sampletimes)
    ix = sortperm(sampletimes; rev=true)
    sampleheights = mst .- sampletimes[ix]
    samplestates = samplestates[ix]
    ushs = sort!(unique(sampleheights))
    deme2shsdict = [countmap(sampleheights[samplestates .== d]) for d in model.demes]

    msol = solveodes(model; solver=solver)
    all_vars = vcat(model.demes, model.nondemes)
    var_index = Dict(v => i for (i, v) in enumerate(all_vars))
    deme_index = Dict(d => i for (i, d) in enumerate(model.demes))
    sym_params = NamedTuple{Tuple(Symbol.(keys(model.parameters)))}(Tuple(values(model.parameters)))
    str_params = model.parameters

    A = Dict(d => 0.0 for d in model.demes)
    events = Event[]
    sizehint!(events, n * 4)

    birthmigrxns = [r for r in model.birthrxn if r.source != r.recipient]
    eventtypes = vcat(
        fill(COALESCENT, length(model.birthrxn)),
        fill(MIGRATION, length(birthmigrxns)),
        fill(MIGRATION, length(model.migrationrxn))
    )
    eventrxns = vcat(model.birthrxn, birthmigrxns, model.migrationrxn)
    nevents = length(eventrxns)
    eventrates = zeros(Float64, nevents)

    function compute_eventrates!(rates::Vector{Float64}, h::Real)
        t_fwd = clamp(mst - h, model.t0, model.tfin)
        u = msol(t_fwd)
        idx = 1
        # 1. Coalescent rates from birth reactions
        for rxn in model.birthrxn
            src = rxn.source
            rec = rxn.recipient
            As = A[src]
            Ar = A[rec]
            Ys = max(As * ytol, u[deme_index[src]], 1e-12)
            Yr = max(Ar * ytol, u[deme_index[rec]], 1e-12)
            f_val = max(0.0, _eval_rxn_rate(rxn, u, var_index, sym_params, str_params, t_fwd))
            if src == rec
                rates[idx] = (As * max(0.0, Ar - 1.0) / (Ys * Yr)) * f_val
            else
                rates[idx] = (As * Ar / (Ys * Yr)) * f_val
            end
            idx += 1
        end
        # 2. Birth-migration rates (cross-deme transmission where parent is not extant)
        for rxn in birthmigrxns
            src = rxn.source
            rec = rxn.recipient
            As = A[src]
            Ar = A[rec]
            u_src = u[deme_index[src]]
            Ys = max(As * ytol, u_src, 1e-12)
            Yr = max(Ar * ytol, u[deme_index[rec]], 1e-12)
            f_val = max(0.0, _eval_rxn_rate(rxn, u, var_index, sym_params, str_params, t_fwd))
            p_not_sampled = u_src > As ? clamp((u_src - As) / max(Ys, 1.0), 0.0, 1.0) : 0.0
            rates[idx] = (Ar / Yr) * p_not_sampled * f_val
            idx += 1
        end
        # 3. Migration rates
        for rxn in model.migrationrxn
            rec = rxn.recipient
            Ar = A[rec]
            Yr = max(Ar * ytol, u[deme_index[rec]], 1e-12)
            g_val = max(0.0, _eval_rxn_rate(rxn, u, var_index, sym_params, str_params, t_fwd))
            rates[idx] = (Ar / Yr) * g_val
            idx += 1
        end
        for i in eachindex(rates)
            if isnan(rates[i]) || rates[i] < 0.0
                rates[i] = 0.0
            end
        end
        return nothing
    end

    # Add initial samples at height 0.0
    ixsampleheight = Ref(1)
    if !isempty(ushs) && ushs[1] == 0.0
        for (ideme, deme) in enumerate(model.demes)
            shsdict = deme2shsdict[ideme]
            if haskey(shsdict, 0.0)
                newsamps = shsdict[0.0]
                A[deme] += newsamps
                for _ in 1:newsamps
                    push!(events, Event(SAMPLE, 0.0, deme, deme))
                end
            end
        end
        ixsampleheight[] = min(length(ushs) + 1, 2)
    end

    cou = Ref(rand(rng))

    function coodes!(du, u, p, h)
        if sum(values(A)) <= 1.0 && ixsampleheight[] > length(ushs)
            du[1] = 0.0
        else
            compute_eventrates!(eventrates, h)
            du[1] = max(0.0, sum(eventrates))
        end
        return nothing
    end

    function sampcondition(u, t, integrator)
        return ixsampleheight[] <= length(ushs) && t == ushs[ixsampleheight[]]
    end

    function sampaffect!(integrator)
        h_curr = ushs[ixsampleheight[]]
        for (ideme, deme) in enumerate(model.demes)
            shsdict = deme2shsdict[ideme]
            if haskey(shsdict, h_curr)
                newsamps = shsdict[h_curr]
                A[deme] += newsamps
                for _ in 1:newsamps
                    push!(events, Event(SAMPLE, h_curr, deme, deme))
                end
            end
        end
        ixsampleheight[] += 1
        return nothing
    end

    function eventcondition(u, t, integrator)
        return exp(-u[1]) - cou[]
    end

    function eventaffect!(integrator)
        compute_eventrates!(eventrates, integrator.t)
        tot_r = sum(eventrates)
        if tot_r > 0.0 && sum(values(A)) > 1.0
            we = sample(rng, 1:nevents, Weights(eventrates))
            rxntype = eventtypes[we]
            rxn = eventrxns[we]
            if rxntype == COALESCENT
                A[rxn.recipient] = max(0.0, A[rxn.recipient] - 1.0)
                push!(events, Event(COALESCENT, integrator.t, rxn.source, rxn.recipient))
            elseif rxntype == MIGRATION
                if A[rxn.recipient] >= 1.0
                    A[rxn.recipient] -= 1.0
                    A[rxn.source] += 1.0
                    push!(events, Event(MIGRATION, integrator.t, rxn.source, rxn.recipient))
                end
            end
        end
        integrator.u[1] = 0.0
        cou[] = rand(rng)
        if sum(values(A)) <= 1.0 && ixsampleheight[] > length(ushs)
            terminate!(integrator)
        end
        return nothing
    end

    sampcb = DiscreteCallback(sampcondition, sampaffect!)
    eventcb = ContinuousCallback(eventcondition, eventaffect!)
    cbs = CallbackSet(sampcb, eventcb)

    pr = ODEProblem(coodes!, [0.0], (0.0, mst - model.t0))
    tstops = [h for h in ushs if h > 0.0]
    solve(pr, solver; callback=cbs, tstops=tstops, maxiters=10^6)

    return SimTree(events, model; computedescendants=computedescendants, rng=rng)
end

"""
    tonewick(o::SimTree) -> String

Convert a `SimTree` to a Newick string (from `Coalescent.jl`).
"""
function tonewick(o::SimTree)
    return write_newick(simtree_to_phylotree(o))
end

"""
    distancematrix(t::SimTree) -> Matrix{Float64}

Compute the pairwise patristic distance matrix between tips of `t::SimTree` (from `Coalescent.jl`).
"""
function distancematrix(t::SimTree)::Matrix{Float64}
    isnothing(t.descendants) && throw(ArgumentError("SimTree must be generated with computedescendants=true."))
    d = zeros(Float64, t.n, t.n)
    dnt = zeros(Float64, t.n + t.nNode, t.n)
    ponodes = sortperm(Float64.(t.heights[(t.n+1):end])) .+ t.n
    for a in ponodes
        u, v, ieu, iev = t.daughters[a - t.n]
        elu = Float64(t.edgelength[ieu])
        elv = Float64(t.edgelength[iev])
        if u > t.n
            descu = findall(t.descendants[u - t.n, :])
            dnt[a, descu] .= elu .+ dnt[u, descu]
        else
            descu = [u]
            dnt[a, u] = elu + dnt[u, u]
        end
        if v > t.n
            descv = findall(t.descendants[v - t.n, :])
            dnt[a, descv] .= elv .+ dnt[v, descv]
        else
            descv = [v]
            dnt[a, v] = elv + dnt[v, v]
        end
        for iu in descu, iv in descv
            dist_uv = dnt[a, iu] + dnt[a, iv]
            d[iu, iv] = dist_uv
            d[iv, iu] = dist_uv
        end
    end
    return d
end

"""
    simtree_to_phylotree(st::SimTree; T_sample::Float64=12.85) -> PhyloTree

Convert a `Coalescent.jl` `SimTree` into a `PhyloTree` (root = `n + 1`, APE node numbering).
"""
function simtree_to_phylotree(st::SimTree; T_sample::Float64=12.85, demes::Union{Nothing,Vector{String}}=nothing)
    n = st.n
    nNode = st.nNode
    nedges = length(st.parent)
    # In SimTree, root is node n + 1 (since a starts at n + nNode and decrements to n + 1)
    edge = zeros(Int, nedges, 2)
    edge_length = zeros(Float64, nedges)
    for e in 1:nedges
        edge[e, 1] = Int(st.parent[e])
        edge[e, 2] = Int(st.child[e])
        edge_length[e] = max(1e-8, Float64(st.edgelength[e]))
    end
    tip_label = [String( something(st.tiplabs[i], "t.$i") ) for i in 1:n]
    deme_idx_fn(lbl::AbstractString, dm::AbstractString) = begin
        if demes !== nothing && !isempty(dm)
            k = findfirst(==(dm), demes)
            k !== nothing && return k
        end
        return infer_tip_state(lbl, dm)
    end
    tip_state = [deme_idx_fn(tip_label[i], st.demes !== nothing ? st.demes[i] : "") for i in 1:n]
    node_state = zeros(Int, n + nNode)
    node_state[1:n] .= tip_state
    if st.demes !== nothing
        for u in (n+1):(n+nNode)
            if isassigned(st.demes, u)
                node_state[u] = deme_idx_fn("", st.demes[u])
            end
        end
    end
    node_time = [T_sample - Float64(something(st.heights[u], 0.0)) for u in 1:(n+nNode)]
    return PhyloTree(
        edge, edge_length, tip_label, tip_state, node_state, node_time,
        nNode, st.true_jumps_12, st.true_jumps_21,
        copy(st.true_jump_times_12), copy(st.true_jump_times_21)
    )
end
