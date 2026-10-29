# ==============================================================================
# phydynr_phyland.jl
# Complete Julia implementation of phydynR (https://github.com/emvolz-phylodynamics/phydynR)
# and phyland (https://github.com/emvolz-phylodynamics/phyland):
# - DatedTree, TFGY, build_demographic_process
# - eventTimes2extant, update_alpha0, rco_finite_size_correction2!
# - solvePikL0 (PL1), solvePikL1 (PL2), solveQAL (QL), solveQfwd0
# - colik / colik_pik / colik_pik_fgy
# - ace (forward-backward Ancestral Character Estimation under structured coalescent)
# - gen_island_model, tips2states, phylandml, confint_phylandml, boot_phylandml
# ==============================================================================

"""
    TFGY

Reverse-time trajectory of demographic matrices `(times, Fs, Gs, Ys, deaths, demes)`
used by `phydynR` (`colik_pik_fgy` and `ace`).
`times` is ordered in decreasing forward time from `t1` (`maxSampleTime`) down to `t0`.
"""
struct TFGY
    times::Vector{Float64}
    Fs::Vector{Matrix{Float64}}
    Gs::Vector{Matrix{Float64}}
    Ys::Vector{Vector{Float64}}
    deaths::Vector{Vector{Float64}}
    demes::Vector{String}
end

"""
    DatedTree

Dated phylogenetic tree with tip sample times and tip deme state probabilities (`phydynR::DatedTree`).
"""
mutable struct DatedTree
    tree::PhyloTree
    n::Int
    m::Int
    Nnode::Int
    demes::Vector{String}
    sampleTimes::Vector{Float64}
    sampleStates::Matrix{Float64}
    maxSampleTime::Float64
    heights::Vector{Float64}
    parentheight::Vector{Float64}
    maxHeight::Float64
    daughters::Matrix{Int}
    parent::Vector{Int}
    lstates::Matrix{Float64}
    mstates::Matrix{Float64}
    ustates::Matrix{Float64}
    acestates::Matrix{Float64}
    coalescentRates::Vector{Float64}
    A::Matrix{Float64}
    lnS::Vector{Float64}
    lnr::Vector{Float64}
    ih::Vector{Float64}
end

function Base.show(io::IO, bdt::DatedTree)
    print(io, "DatedTree(n=$(bdt.n) tips, m=$(bdt.m) demes $(bdt.demes), maxSampleTime=$(round(bdt.maxSampleTime, digits=3)), maxHeight=$(round(bdt.maxHeight, digits=3)))")
end

"""
    DatedTree(tree::PhyloTree; sampleTimes=nothing, sampleStates=nothing, demes=["Y1", "Y2"], minEdgeLength=1e-6) -> DatedTree

Construct a `phydynR`-compatible `DatedTree` from a `PhyloTree`.
"""
function DatedTree(
    tree_in::PhyloTree;
    sampleTimes::Union{Nothing,AbstractVector{<:Real}} = nothing,
    sampleStates::Union{Nothing,AbstractMatrix{<:Real}} = nothing,
    demes::Vector{String} = ["Y1", "Y2"],
    minEdgeLength::Float64 = 1e-6
)
    bt = collapse_singletons(tree_in)
    bt.edge_length .= max.(bt.edge_length, minEdgeLength)
    n = ntips(bt)
    Nnode = bt.nnode
    ntot = n + Nnode

    dists_from_root = node_heights(bt)
    if sampleTimes === nothing
        if !isempty(bt.node_time) && all(!isnan, bt.node_time[1:n]) && maximum(bt.node_time[1:n]) > 0.0
            sts = Float64.(bt.node_time[1:n])
        else
            sts = Float64.(dists_from_root[1:n])
        end
    else
        sts = Float64.(sampleTimes)
    end

    if sampleStates === nothing
        m = length(demes)
        ssts = zeros(Float64, n, m)
        for i in 1:n
            st = clamp(bt.tip_state[i], 1, m)
            ssts[i, st] = 1.0
        end
    else
        ssts = Float64.(sampleStates)
        m = size(ssts, 2)
        if length(demes) != m
            demes = ["deme_$k" for k in 1:m]
        end
    end

    max_st = maximum(sts)
    max_root_dist = maximum(dists_from_root[1:n])
    heights = [max(0.0, max_root_dist - dists_from_root[u]) for u in 1:ntot]
    # Round tiny tip heights for homochronous tips
    for i in 1:n
        if abs(heights[i]) < 1e-5
            heights[i] = 0.0
        end
    end
    max_h = maximum(heights)

    daughters = fill(0, ntot, 2)
    parent = fill(0, ntot)
    parentheight = fill(NaN, ntot)
    child_count = zeros(Int, ntot)
    for e in 1:size(bt.edge, 1)
        p = bt.edge[e, 1]
        c = bt.edge[e, 2]
        parent[c] = p
        parentheight[c] = heights[p]
        child_count[p] += 1
        if child_count[p] <= 2
            daughters[p, child_count[p]] = c
        end
    end

    lstates = fill(NaN, m, ntot)
    lstates[:, 1:n] .= transpose(ssts)
    mstates = copy(lstates)
    ustates = zeros(Float64, m, ntot)
    acestates = fill(NaN, m, ntot)
    coalescentRates = zeros(Float64, ntot)

    return DatedTree(
        bt, n, m, Nnode, demes, sts, ssts, max_st,
        heights, parentheight, max_h, daughters, parent,
        lstates, mstates, ustates, acestates, coalescentRates,
        zeros(Float64, 0, m), Float64[], Float64[], Float64[]
    )
end

DatedTree(st::SimTree; kwargs...) = DatedTree(simtree_to_phylotree(st); kwargs...)

"""
    eventTimes2extant(eventTimes, nodeheights, parentheights) -> (extant, nodesAtHeight)

Identify extant lineages and nodes at each event height (ported from `phydynR/src/colikModular0.cpp`).
"""
function eventTimes2extant(
    eventTimes::AbstractVector{Float64},
    nodeheights::AbstractVector{Float64},
    parentheights::AbstractVector{Float64}
)
    nev = length(eventTimes)
    nnodes = length(nodeheights)
    extant = [Int[] for _ in 1:nev]
    nodesAtHeight = [Int[] for _ in 1:nev]
    for i in 1:nev
        h = eventTimes[i]
        ext_i = extant[i]
        nah_i = nodesAtHeight[i]
        for u in 1:nnodes
            nh = nodeheights[u]
            ph = parentheights[u]
            if nh == h
                push!(nah_i, u)
            end
            if !isnan(ph)
                if nh <= h && ph > h
                    push!(ext_i, u)
                end
            else
                if nh <= h
                    push!(ext_i, u)
                end
            end
        end
    end
    return extant, nodesAtHeight
end

"""
    update_alpha0(pu, pv, F, Y, A; min_y=1.0) -> (pa, corate)

Compute ancestral node state distribution `pa` and instantaneous coalescence rate `corate`
at a coalescent node with daughters `(u, v)` (ported from `phydynR/src/colikModular0.cpp`).
"""
function update_alpha0(
    pu::AbstractVector{Float64},
    pv::AbstractVector{Float64},
    F::AbstractMatrix{Float64},
    Y::AbstractVector{Float64},
    A::AbstractVector{Float64};
    min_y::Float64 = 1e-6
)
    m = length(Y)
    Y_eff = [max(Y[k], A[k], min_y) for k in 1:m]
    pu_Y = pu ./ Y_eff
    pv_Y = pv ./ Y_eff
    pa = pu_Y .* (F * pv_Y) .+ pv_Y .* (F * pu_Y)
    corate = sum(pa)
    if corate > 0.0 && isfinite(corate)
        pa ./= corate
    else
        fill!(pa, 1.0 / m)
        corate = 0.0
    end
    return pa, corate
end

"""
    rco_finite_size_correction2!(alpha, p_a, A, extantLines, mstates)

Apply Volz (2012) finite-size correction to non-coalescing extant lineages
(ported from `phydynR/src/colikModular0.cpp`).
"""
function rco_finite_size_correction2!(
    alpha::Int,
    p_a::AbstractVector{Float64},
    A::AbstractVector{Float64},
    extantLines::AbstractVector{Int},
    mstates::AbstractMatrix{Float64}
)
    m = size(mstates, 1)
    for u in extantLines
        if u != alpha
            p_u = view(mstates, :, u)
            denom = [max(A[k] - p_u[k], 1e-12) for k in 1:m]
            rterm = p_a ./ denom
            rho = A ./ denom
            lterm = dot(rho, p_a)
            for k in 1:m
                p_u[k] = p_u[k] * max(lterm - rterm[k], 0.0)
            end
            s = sum(p_u)
            if s > 0.0
                p_u ./= s
            end
        end
    end
    return mstates
end

@inline function _fgy_index_at_height(h::Float64, hres::Int, treeT::Float64)
    idx = floor(Int, hres * h / max(treeT, 1e-12)) + 1
    return clamp(idx, 1, hres)
end

"""
    solvePikL0(tfgy::TFGY, h0::Float64, h1::Float64, pik0::AbstractMatrix{Float64}, step_size_res::Real=10)

Integrate lineage state probabilities `Pik` and cumulative hazard `L` over `[h0, h1]`
under the `PL1` likelihood approximation (ported from `phydynR/src/dPikL0.cpp`).
"""
function solvePikL0(
    tfgy::TFGY,
    h0::Float64,
    h1::Float64,
    pik0::AbstractMatrix{Float64},
    step_size_res::Real = 10;
    min_y::Float64 = 1.0
)
    m, nextant = size(pik0)
    if h1 <= h0 || nextant == 0
        return (pik = copy(pik0), L = 0.0)
    end
    hres = length(tfgy.times)
    treeT = abs(tfgy.times[1] - tfgy.times[end])
    dt_grid = hres > 1 ? abs(tfgy.times[2] - tfgy.times[1]) : (h1 - h0)
    nsteps = max(2, ceil(Int, (h1 - h0) / min(dt_grid / max(1.0, Float64(step_size_res)), (h1 - h0) / max(1.0, Float64(step_size_res)))))
    dh = (h1 - h0) / nsteps

    Pik = copy(pik0)
    dPik = zeros(Float64, m, nextant)
    R = zeros(Float64, m, m)
    A = zeros(Float64, m)
    L = 0.0

    for step in 1:nsteps
        h_mid = h0 + (step - 0.5) * dh
        idx = _fgy_index_at_height(h_mid, hres, treeT)
        F = tfgy.Fs[idx]
        G = tfgy.Gs[idx]
        Y = tfgy.Ys[idx]
        @inbounds for k in 1:m
            A[k] = 0.0
            for z in 1:nextant
                A[k] += Pik[k, z]
            end
        end
        @inbounds for l in 1:m
            Yl = max(Y[l], min_y)
            for k in 1:m
                R[k, l] = k == l ? 0.0 : (F[k, l] + G[k, l]) / Yl
            end
        end
        @inbounds for k in 1:m
            s = 0.0
            for l in 1:m
                l != k && (s += R[l, k])
            end
            R[k, k] = -s
        end

        mul!(dPik, R, Pik)
        dL = 0.0
        @inbounds for z in 1:nextant
            for k in 1:m
                Yk = max(Y[k], min_y)
                pkz_Yk = Pik[k, z] / Yk
                for l in 1:m
                    Yl = max(Y[l], min_y)
                    amp_Yl = clamp(max(0.0, A[l] - Pik[l, z]) / Yl, 0.0, 1e6)
                    term = pkz_Yk * F[k, l] * amp_Yl
                    dPik[k, z] -= term
                    dL += term * 0.5
                end
            end
        end
        Pik .= max.(0.0, Pik .+ dh .* dPik)
        L += dh * dL
    end
    return (pik = Pik, L = L)
end

"""
    solvePikL1(tfgy::TFGY, h0::Float64, h1::Float64, pik0::AbstractMatrix{Float64}, step_size_res::Real=10)

Integrate lineage state probabilities `Pik` and cumulative hazard `L` over `[h0, h1]`
under the `PL2` likelihood approximation (ported from `phydynR/src/dPikL1.cpp`).
"""
function solvePikL1(
    tfgy::TFGY,
    h0::Float64,
    h1::Float64,
    pik0::AbstractMatrix{Float64},
    step_size_res::Real = 10;
    min_y::Float64 = 1e-6
)
    m, nextant = size(pik0)
    if h1 <= h0 || nextant == 0
        return (pik = copy(pik0), L = 0.0)
    end
    hres = length(tfgy.times)
    treeT = abs(tfgy.times[1] - tfgy.times[end])
    dt_grid = hres > 1 ? abs(tfgy.times[2] - tfgy.times[1]) : (h1 - h0)
    nsteps = max(2, ceil(Int, (h1 - h0) / min(dt_grid / max(1.0, Float64(step_size_res)), (h1 - h0) / max(1.0, Float64(step_size_res)))))
    dh = (h1 - h0) / nsteps

    Pik = copy(pik0)
    nPik = zeros(Float64, m, nextant)
    dPik = zeros(Float64, m, nextant)
    R = zeros(Float64, m, m)
    phi = zeros(Float64, m, m)
    nA = zeros(Float64, m)
    u_vec = zeros(Float64, m)
    L = 0.0

    for step in 1:nsteps
        h_mid = h0 + (step - 0.5) * dh
        idx = _fgy_index_at_height(h_mid, hres, treeT)
        F = tfgy.Fs[idx]
        G = tfgy.Gs[idx]
        Y = tfgy.Ys[idx]

        fill!(nA, 0.0)
        @inbounds for z in 1:nextant
            s = 0.0
            for k in 1:m
                s += max(0.0, Pik[k, z])
            end
            if s > 0.0
                for k in 1:m
                    v = max(0.0, Pik[k, z]) / s
                    nPik[k, z] = v
                    nA[k] += v
                end
            else
                for k in 1:m
                    nPik[k, z] = 0.0
                end
            end
        end

        @inbounds for k in 1:m
            Yk = max(Y[k], min_y)
            u_vec[k] = Y[k] > nA[k] ? clamp((Y[k] - nA[k]) / max(Yk, 1.0), 0.0, 1.0) : 0.0
        end

        @inbounds for l in 1:m
            Yl = max(Y[l], nA[l], min_y)
            for k in 1:m
                Yk = max(Y[k], nA[k], min_y)
                phi[k, l] = F[k, l] / (Yk * Yl)
                if k == l
                    R[k, l] = 0.0
                else
                    raw_r = (F[k, l] * u_vec[k] + G[k, l]) / Yl
                    cap_r = nA[l] > 1e-12 ? max(0.0, Y[k] - nA[k]) / (dh * nA[l]) : raw_r
                    R[k, l] = min(raw_r, cap_r)
                end
            end
        end
        @inbounds for k in 1:m
            s = 0.0
            for l in 1:m
                l != k && (s += R[l, k])
            end
            R[k, k] = -s
        end

        # Exact matrix exponential update exp(R * dh) * nPik (unconditionally stable & probability-conserving)
        P_step = exp(dh .* R)
        # Midpoint state for hazard integration
        P_half = exp((0.5 * dh) .* R)
        mul!(dPik, P_half, nPik)

        fill!(nA, 0.0)
        @inbounds for z in 1:nextant
            for k in 1:m
                nA[k] += dPik[k, z]
            end
        end

        dL = 0.0
        @inbounds for z in 1:nextant
            s_col = sum(view(Pik, :, z))
            s_col <= 0.0 && continue
            for k in 1:m
                for l in 1:m
                    phi_sym = phi[k, l] + phi[l, k]
                    rem_l = max(0.0, nA[l] - dPik[l, z])
                    dL += 0.5 * dPik[k, z] * phi_sym * rem_l
                end
            end
        end

        mul!(Pik, P_step, nPik)
        Pik .= max.(0.0, Pik)
        L += dh * dL
    end
    return (pik = Pik, L = L)
end

"""
    solveQAL(tfgy::TFGY, h0::Float64, h1::Float64, A0::AbstractVector{Float64}, L0::Float64=0.0; step_size_res::Real=10)

Solve transition matrix `Q(h0, h1)`, aggregate lineage counts `A(h1)`, and cumulative hazard `L`
under the `QL` approximation (ported from `phydynR/src/solveQALboost0.cpp` and fixed `dQAL-6.3.c`).
"""
function solveQAL(
    tfgy::TFGY,
    h0::Float64,
    h1::Float64,
    A0::AbstractVector{Float64},
    L0::Float64 = 0.0;
    step_size_res::Real = 10,
    min_y::Float64 = 1e-6
)
    m = length(A0)
    Q = Matrix{Float64}(I, m, m)
    A = copy(A0)
    Atotal = sum(A0)
    L = L0
    if h1 <= h0 || Atotal <= 0.0
        return (Q = Q, A = A, L = L)
    end
    hres = length(tfgy.times)
    treeT = abs(tfgy.times[1] - tfgy.times[end])
    dt_grid = hres > 1 ? abs(tfgy.times[2] - tfgy.times[1]) : (h1 - h0)
    nsteps = max(2, ceil(Int, (h1 - h0) / min(dt_grid / max(1.0, Float64(step_size_res)), (h1 - h0) / max(1.0, Float64(step_size_res)))))
    dh = (h1 - h0) / nsteps

    dQ = zeros(Float64, m, m)
    dA = zeros(Float64, m)
    pair_factor = Atotal > 1.0 ? ((Atotal - 1.0) / Atotal) : 0.0
    for step in 1:nsteps
        h_mid = h0 + (step - 0.5) * dh
        idx = _fgy_index_at_height(h_mid, hres, treeT)
        F = tfgy.Fs[idx]
        G = tfgy.Gs[idx]
        Y = tfgy.Ys[idx]

        sumA = sum(A)
        r = sumA > 1e-12 ? (Atotal / sumA) : 1.0
        fill!(dA, 0.0)
        fill!(dQ, 0.0)
        dL = 0.0
        @inbounds for k in 1:m
            Yk = max(Y[k], min_y)
            ak = clamp(r * A[k] / Yk, 0.0, 1.0)
            for l in 1:m
                Yl = max(Y[l], min_y)
                al = clamp(r * A[l] / Yl, 0.0, 1.0)
                if l != k
                    dA[k] += ((1.0 - ak) * F[k, l] + G[k, l]) * al - ((1.0 - al) * F[l, k] + G[l, k]) * ak
                    for z in 1:m
                        dQ[k, z] += (F[k, l] * (1.0 - ak) + G[k, l]) * (Q[l, z] / Yl) - (F[l, k] * (1.0 - al) + G[l, k]) * (Q[k, z] / Yk)
                    end
                    dL += max(0.0, ak * al * F[k, l])
                else
                    dL += max(0.0, pair_factor * ((r * A[k]) / Yk)^2 * F[k, k])
                end
            end
        end
        A .= max.(0.0, A .+ dh .* dA)
        Q .= max.(0.0, Q .+ dh .* dQ)
        for z in 1:m
            cs = sum(view(Q, :, z))
            if cs > 0.0
                Q[:, z] ./= cs
            end
        end
        L += dh * dL
    end
    return (Q = Q, A = A, L = L)
end

"""
    solveQfwd0(tfgy::TFGY, m::Int, h1::Float64, h0::Float64; min_y::Float64=1.0) -> Matrix{Float64}

Solve forward-in-time transition matrix `Q` from height `h1` down to height `h0` (`h1 > h0`)
for ancestral character estimation `ace` (ported from `phydynR/src/dQfwd.cpp`).
"""
function solveQfwd0(
    tfgy::TFGY,
    m::Int,
    h1::Float64,
    h0::Float64;
    min_y::Float64 = 1.0
)::Matrix{Float64}
    Q = Matrix{Float64}(I, m, m)
    h1 <= h0 && return Q
    hres = length(tfgy.times)
    treeT = abs(tfgy.times[1] - tfgy.times[end])
    t0 = treeT - h1
    t1 = treeT - h0
    nsteps = 100
    dt = (t1 - t0) / nsteps
    dQ = zeros(Float64, m, m)

    for step in 1:nsteps
        t_mid = t0 + (step - 0.5) * dt
        h_mid = treeT - t_mid
        idx = _fgy_index_at_height(h_mid, hres, treeT)
        F = tfgy.Fs[idx]
        G = tfgy.Gs[idx]
        Y = tfgy.Ys[idx]
        death = tfgy.deaths[idx]
        fill!(dQ, 0.0)
        @inbounds for z in 1:m
            for k in 1:m
                for l in 1:m
                    if k != l
                        if Q[l, z] > 0.0
                            dQ[k, z] += (0.5 * F[l, k] + G[l, k]) * Q[l, z] / max(Q[l, z], Y[l], min_y)
                        end
                        if Q[k, z] > 0.0
                            dQ[k, z] -= (0.5 * F[k, l] + G[k, l]) * Q[k, z] / max(Q[k, z], Y[k], min_y)
                        end
                    end
                end
                if Q[k, z] > 0.0
                    dQ[k, z] -= Q[k, z] * death[k] / max(Q[k, z], Y[k], min_y)
                end
            end
        end
        Q .= max.(0.0, Q .+ dt .* dQ)
    end
    for z in 1:m
        cs = sum(view(Q, :, z))
        if cs > 0.0
            Q[:, z] ./= cs
        end
    end
    return Q
end

# Order concurrent internal nodes so daughters are processed before parents
function _new_nodes_order(nn::Vector{Int}, extantLines::Vector{Int}, bdt::DatedTree)
    length(nn) <= 1 && return nn
    o = Int[]
    nrep = 0
    avail = union(collect(1:bdt.n), extantLines)
    while length(o) < length(nn)
        for a in setdiff(nn, o)
            d1 = bdt.daughters[a, 1]
            d2 = bdt.daughters[a, 2]
            if (d1 in avail || d1 in o) && (d2 in avail || d2 in o)
                push!(o, a)
            end
        end
        nrep += 1
        if nrep > 1 + length(nn)
            o = copy(nn)
            break
        end
    end
    return o
end

"""
    colik_pik_fgy(bdt_in::DatedTree, tfgy::TFGY; kwargs...)

Compute the structured coalescent log-likelihood given a `DatedTree` and a `TFGY` demographic history
(ported from `phydynR/R/colik.pik0.R` and `colikModular.R`).
Supports `likelihood = :PL2` (default), `:PL1`, and `:QL`.
"""
function colik_pik_fgy(
    bdt_in::DatedTree,
    tfgy::TFGY;
    timeOfOriginBoundaryCondition::Bool = true,
    maxHeight::Float64 = Inf,
    forgiveAgtY::Float64 = 1.0,
    AgtY_penalty::Float64 = 0.0,
    returnTree::Bool = false,
    step_size_res::Real = 10,
    likelihood::Symbol = :PL2,
    finiteSizeCorrection::Bool = false
)
    bdt = returnTree ? deepcopy(bdt_in) : bdt_in
    t0 = tfgy.times[end]
    if min(bdt.maxHeight, maxHeight) > (bdt.maxSampleTime - t0) && timeOfOriginBoundaryCondition
        return returnTree ? (loglik = -Inf, tree = bdt) : -Inf
    end

    hres = length(tfgy.times)
    treeT = abs(tfgy.times[1] - tfgy.times[end])
    get_fgy(h) = begin
        ih = _fgy_index_at_height(h, hres, treeT)
        (F = tfgy.Fs[ih], G = tfgy.Gs[ih], Y = tfgy.Ys[ih])
    end

    ntot = bdt.n + bdt.Nnode
    fill!(bdt.lstates, NaN)
    bdt.lstates[:, 1:bdt.n] .= transpose(bdt.sampleStates)
    bdt.mstates .= bdt.lstates
    fill!(bdt.ustates, 0.0)

    eventTimes = sort!(unique(bdt.heights))
    if isfinite(maxHeight)
        filter!(h -> h <= maxHeight, eventTimes)
    end
    extantAtEvent_list, nodesAtHeight = eventTimes2extant(eventTimes, bdt.heights, bdt.parentheight)

    A_trace = Vector{Float64}[]
    empty!(bdt.lnS)
    empty!(bdt.lnr)
    empty!(bdt.ih)

    L = 0.0
    loglik = 0.0

    for ih in 1:(length(eventTimes) - 1)
        h0 = eventTimes[ih]
        h1 = eventTimes[ih + 1]
        fgy = get_fgy(h1)
        extantLines = extantAtEvent_list[ih]
        isempty(extantLines) && continue

        if likelihood === :QL
            A0 = vec(sum(view(bdt.mstates, :, extantLines); dims=2))
            qal = solveQAL(tfgy, h0, h1, A0, 0.0; step_size_res=step_size_res)
            Q = qal.Q
            L += qal.L
            for u in extantLines
                col = Q * view(bdt.mstates, :, u)
                col .= clamp.(col, 0.0, 1.0)
                s = sum(col)
                if s > 0.0
                    bdt.mstates[:, u] .= col ./ s
                end
            end
        else
            # PL2 or PL1 (append dummy column of zeros like phydynR's cbind(mstates[,extantLines], 0))
            pik0 = hcat(view(bdt.mstates, :, extantLines), zeros(Float64, bdt.m))
            pikl = likelihood === :PL2 ?
                solvePikL1(tfgy, h0, h1, pik0, step_size_res) :
                solvePikL0(tfgy, h0, h1, pik0, step_size_res)
            pik1 = max.(pikl.pik[:, 1:length(extantLines)], 0.0)
            for j in 1:length(extantLines)
                s = sum(view(pik1, :, j))
                if s > 0.0
                    pik1[:, j] ./= s
                else
                    pik1[:, j] .= 1.0 / bdt.m
                end
            end
            bdt.mstates[:, extantLines] .= pik1
            L += pikl.L
        end

        newNodes = [u for u in nodesAtHeight[ih + 1] if u > bdt.n]
        A = vec(sum(view(bdt.mstates, :, extantLines); dims=2))
        YmA = sum(fgy.Y) - length(extantLines)
        if YmA < 0.0
            if (length(extantLines) / bdt.n) > forgiveAgtY
                L = Inf
            else
                L += L * abs(YmA) * AgtY_penalty
            end
        end

        for alpha in _new_nodes_order(newNodes, extantLines, bdt)
            u = bdt.daughters[alpha, 1]
            v = bdt.daughters[alpha, 2]
            pa, corate = update_alpha0(view(bdt.mstates, :, u), view(bdt.mstates, :, v), fgy.F, fgy.Y, A)
            bdt.coalescentRates[alpha] = corate
            if isnan(L)
                L = Inf
            end
            bdt.lstates[:, alpha] .= pa
            bdt.mstates[:, alpha] .= pa

            push!(A_trace, copy(A))
            push!(bdt.lnS, -L)
            push!(bdt.lnr, corate > 0.0 ? log(corate) : -Inf)
            push!(bdt.ih, h1)

            if corate <= 0.0 || !isfinite(corate) || !isfinite(L)
                return returnTree ? (loglik = -Inf, tree = bdt) : -Inf
            end
            loglik += log(corate) - L
            if !isfinite(loglik)
                return returnTree ? (loglik = -Inf, tree = bdt) : -Inf
            end
            if finiteSizeCorrection
                rco_finite_size_correction2!(alpha, pa, A, extantLines, bdt.mstates)
            end
        end

        if !isempty(newNodes)
            L = 0.0
        end
    end

    if !isempty(A_trace)
        bdt.A = reduce(vcat, [transpose(a) for a in A_trace])
    end

    return returnTree ? (loglik = loglik, tree = bdt) : loglik
end

"""
    build_demographic_process(model::ModelFGY) -> Function

Build a `phydynR`-style demographic process callable `dm(theta, x0, t0, t1; res=1000)`
from a `Coalescent.jl` `ModelFGY`.
"""
function build_demographic_process(base_model::ModelFGY)
    function dm(
        theta::AbstractDict = Dict{String,Float64}(),
        x0::Union{Nothing,AbstractDict} = nothing,
        t0::Real = base_model.t0,
        t1::Real = base_model.tfin;
        res::Int = 1000,
        solver = Tsit5()
    )::TFGY
        mod = deepcopy(base_model)
        mod.t0 = Float64(t0)
        mod.tfin = Float64(t1)
        for (k, v) in theta
            mod.parameters[String(k)] = Float64(v)
        end
        if x0 !== nothing
            for (k, v) in x0
                mod.initial[String(k)] = Float64(v)
            end
        end

        sol = solveodes(mod; solver=solver, res=res)
        all_vars = vcat(mod.demes, mod.nondemes)
        var_index = Dict(v => i for (i, v) in enumerate(all_vars))
        deme_index = Dict(d => i for (i, d) in enumerate(mod.demes))
        sym_params = NamedTuple{Tuple(Symbol.(keys(mod.parameters)))}(Tuple(values(mod.parameters)))
        str_params = mod.parameters
        m = mod.numberdemes

        times_rev = reverse(sol.t)
        nt = length(times_rev)
        Fs = Vector{Matrix{Float64}}(undef, nt)
        Gs = Vector{Matrix{Float64}}(undef, nt)
        Ys = Vector{Vector{Float64}}(undef, nt)
        deaths = Vector{Vector{Float64}}(undef, nt)

        for (idx, t_val) in enumerate(times_rev)
            u = sol.u[nt - idx + 1]
            F = zeros(Float64, m, m)
            G = zeros(Float64, m, m)
            Y = [max(0.0, u[i]) for i in 1:m]
            D = zeros(Float64, m)
            for rxn in mod.birthrxn
                r = max(0.0, _eval_rxn_rate(rxn, u, var_index, sym_params, str_params, t_val))
                F[deme_index[rxn.source], deme_index[rxn.recipient]] += r
            end
            for rxn in mod.migrationrxn
                r = max(0.0, _eval_rxn_rate(rxn, u, var_index, sym_params, str_params, t_val))
                G[deme_index[rxn.source], deme_index[rxn.recipient]] += r
            end
            for rxn in mod.deathrxn
                r = max(0.0, _eval_rxn_rate(rxn, u, var_index, sym_params, str_params, t_val))
                D[deme_index[rxn.source]] += r
            end
            Fs[idx] = F
            Gs[idx] = G
            Ys[idx] = Y
            deaths[idx] = D
        end

        return TFGY(times_rev, Fs, Gs, Ys, deaths, copy(mod.demes))
    end
    return dm
end

"""
    colik(bdt::DatedTree, theta, dm::Function; x0=nothing, t0=0.0, res=1000, likelihood=:PL2, kwargs...)

Compute the `phydynR` structured coalescent log-likelihood (alias `colik_pik`).
"""
function colik(
    bdt::DatedTree,
    theta::AbstractDict,
    dm::Function;
    x0::Union{Nothing,AbstractDict} = nothing,
    t0::Real = 0.0,
    res::Int = 1000,
    timeOfOriginBoundaryCondition::Bool = true,
    maxHeight::Float64 = Inf,
    forgiveAgtY::Float64 = 1.0,
    AgtY_penalty::Float64 = 0.0,
    returnTree::Bool = false,
    step_size_res::Real = 10,
    likelihood::Symbol = :PL2,
    finiteSizeCorrection::Bool = false
)
    tfgy = dm(theta, x0, Float64(t0), bdt.maxSampleTime; res=res)
    return colik_pik_fgy(
        bdt, tfgy;
        timeOfOriginBoundaryCondition=timeOfOriginBoundaryCondition,
        maxHeight=maxHeight,
        forgiveAgtY=forgiveAgtY,
        AgtY_penalty=AgtY_penalty,
        returnTree=returnTree,
        step_size_res=step_size_res,
        likelihood=likelihood,
        finiteSizeCorrection=finiteSizeCorrection
    )
end

"""Alias for [`colik`](@ref), the demographic-process structured coalescent log-likelihood."""
const colik_pik = colik

"""
    ace(bdt::DatedTree, theta::AbstractDict, dm::Function; x0=nothing, t0=0.0, res=1000, likelihood=:PL2, kwargs...) -> DatedTree

Forward-backward Ancestral Character Estimation (`phydynR::ace`) under a structured
coalescent demographic model. Populates `bdt.acestates` (`m × (n + Nnode)`) with posterior
ancestral state probabilities at all nodes.
"""
function ace(
    bdt_in::DatedTree,
    theta::AbstractDict,
    dm::Function;
    x0::Union{Nothing,AbstractDict} = nothing,
    t0::Real = 0.0,
    res::Int = 1000,
    timeOfOriginBoundaryCondition::Bool = true,
    maxHeight::Float64 = Inf,
    forgiveAgtY::Float64 = 1.0,
    AgtY_penalty::Float64 = 1.0,
    step_size_res::Real = 10,
    likelihood::Symbol = :PL2
)::DatedTree
    tfgy = dm(theta, x0, Float64(t0), bdt_in.maxSampleTime; res=res)
    m = bdt_in.m
    res_colik = colik_pik_fgy(
        bdt_in, tfgy;
        timeOfOriginBoundaryCondition=timeOfOriginBoundaryCondition,
        maxHeight=maxHeight,
        forgiveAgtY=forgiveAgtY,
        AgtY_penalty=AgtY_penalty,
        returnTree=true,
        step_size_res=step_size_res,
        likelihood=likelihood
    )
    bdt = res_colik.tree
    eventTimes = sort!(unique(bdt.heights))
    if isfinite(maxHeight)
        filter!(h -> h <= maxHeight, eventTimes)
    end
    extantAtEvent_list, nodesAtHeight = eventTimes2extant(eventTimes, bdt.heights, bdt.parentheight)

    hres = length(tfgy.times)
    treeT = abs(tfgy.times[1] - tfgy.times[end])
    get_fgy(h) = begin
        ih = _fgy_index_at_height(h, hres, treeT)
        (F = tfgy.Fs[ih], G = tfgy.Gs[ih], Y = tfgy.Ys[ih])
    end

    lstates = bdt.lstates
    acestates = fill(NaN, size(lstates))
    root = bdt.n + 1
    acestates[:, root] .= lstates[:, root]
    fgy_root = get_fgy(maximum(bdt.heights))
    x = transpose(fgy_root.F) * acestates[:, root]
    sx = sum(x)
    x = sx > 0.0 ? (x ./ sx) : fill(1.0 / m, m)
    x = 0.5 .* (x .+ acestates[:, root])
    for d_idx in 1:2
        c = bdt.daughters[root, d_idx]
        if c > 0
            acestates[:, c] .= x
        end
    end

    for ih in length(eventTimes):-1:2
        h1 = eventTimes[ih]
        h0 = eventTimes[ih - 1]
        extantLines = extantAtEvent_list[ih - 1]
        nodes_h0 = nodesAtHeight[ih - 1]

        if h1 > h0 && !isempty(extantLines)
            QQ = solveQfwd0(tfgy, m, h1, h0)
            acestates[:, extantLines] .= QQ * view(acestates, :, extantLines)
        end

        fgy = get_fgy(h0)
        for u in nodes_h0
            xu = acestates[:, u] .* lstates[:, u]
            sxu = sum(xu)
            if sxu > 0.0 && isfinite(sxu)
                acestates[:, u] .= xu ./ sxu
            elseif !any(isnan, view(lstates, :, u))
                acestates[:, u] .= lstates[:, u]
            end
            v0 = bdt.daughters[u, 1]
            v1 = bdt.daughters[u, 2]
            xd = transpose(fgy.F) * acestates[:, u]
            sxd = sum(xd)
            xd = sxd > 0.0 ? (xd ./ sxd) : copy(acestates[:, u])
            xd = 0.5 .* (xd .+ acestates[:, u])
            v0 > 0 && (acestates[:, v0] .= xd)
            v1 > 0 && (acestates[:, v1] .= xd)
        end
    end

    bdt.acestates = acestates
    return bdt
end

# ==============================================================================
# phyland: Equilibrium Structured Island Coalescent ML & Ancestral Reconstruction
# (ported from https://github.com/emvolz-phylodynamics/phyland/blob/master/R/pgy.R)
# ==============================================================================

"""
    gen_island_model(demes::Vector{String}; xF::Float64=1e4) -> Function

Generate an equilibrium structured island demographic process `dm(theta, x0, t0, t1; res=10)`
(ported from `phyland::.gen.island.model`).
`theta` contains `"Ne1".."Nem"` and `"mu_k_l"` (`k != l`).
"""
function gen_island_model(demes::Vector{String}; xF::Float64 = 1e4)
    m = length(demes)
    function dm_island(
        theta::AbstractDict,
        x0 = nothing,
        t0::Real = 0.0,
        t1::Real = 1.0;
        res::Int = 10
    )::TFGY
        times = collect(reverse(range(Float64(t0), Float64(t1); length=res)))
        b = zeros(Float64, m, m)
        y = zeros(Float64, m)
        for k in 1:m
            nek = max(1e-8, Float64(theta["Ne$k"]))
            b[k, k] = nek * xF * xF / 2.0
            y[k] = nek * xF
        end
        g = zeros(Float64, m, m)
        for k in 1:m, l in 1:m
            if k != l
                mukl = Float64(get(theta, "mu$(l)$(k)", 0.0))
                g[k, l] = y[l] * max(0.0, mukl)
            end
        end
        Fs = [copy(b) for _ in 1:res]
        Gs = [copy(g) for _ in 1:res]
        Ys = [copy(y) for _ in 1:res]
        deaths = [zeros(Float64, m) for _ in 1:res]
        return TFGY(times, Fs, Gs, Ys, deaths, copy(demes))
    end
    return dm_island
end

"""
    tips2states(tips::Vector{String}; delimiter::String="_", index::Union{Nothing,Int}=nothing, regex::Union{Nothing,Regex}=nothing)

Extract deme names and sample-state indicator matrix `ssts` (`n × m`) from tip labels
(ported from `phyland::.tips2states`).
"""
function tips2states(
    tips::Vector{String};
    delimiter::String = "_",
    index::Union{Nothing,Int} = nothing,
    regex::Union{Nothing,Regex} = nothing
)
    tipstates = String[]
    for t in tips
        if regex !== nothing
            m = match(regex, t)
            m === nothing && error("Regex did not match tip label: $t")
            push!(tipstates, String(m.match))
        else
            parts = split(t, delimiter)
            idx = index === nothing ? length(parts) : index
            push!(tipstates, String(parts[clamp(idx, 1, length(parts))]))
        end
    end
    demes = sort!(unique(tipstates))
    m = length(demes)
    ssts = zeros(Float64, length(tips), m)
    for (i, ts) in enumerate(tipstates)
        k = findfirst(==(ts), demes)
        ssts[i, k] = 1.0
    end
    return (ssts = ssts, demes = demes, m = m, tipstates = tipstates)
end

"""
    PhylandFit

Island-model inference result from [`phylandml`](@ref). Contains natural-scale
`coef`, fitted `logcoef` and `param_names`, `loglik`, `demes`, migration `design`,
an ancestral-state table `ace`, the reconstructed `bdt`, and `expected_jumps`
(sums of endpoint posterior products, not stochastic-mapping jump counts).
"""
struct PhylandFit
    coef::Dict{String,Float64}
    logcoef::Vector{Float64}
    param_names::Vector{String}
    loglik::Float64
    demes::Vector{String}
    design::Matrix{Int}
    ace::DataFrame
    bdt::DatedTree
    expected_jumps::Matrix{Float64}
end

"""
    PhylandProfile

Result from [`confint_phylandml`](@ref): parameter name `param`, natural-scale
interval `ci`, log-parameter `grid`, log-likelihoods `ll`, and table `sumdat`.
"""
struct PhylandProfile
    param::String
    ci::Tuple{Float64,Float64}
    grid::Vector{Float64}
    ll::Vector{Float64}
    sumdat::DataFrame
end

"""
    phylandml(tree::PhyloTree; delimiter="_", index=nothing, regex=nothing, design=nothing, start_migrate=NaN, start_Ne=NaN, minBL=1e-6, likelihood=:PL2) -> PhylandFit

Maximum likelihood inference of effective population sizes `Ne_k` and migration rates `mu_{l->k}`
between demes, plus ancestral character estimation (`ace`) (ported from `phyland::phylandml`).
In `design[k, l]`, `0` fixes a migration rate to `0.0`, while positive integers group shared migration rates.
"""
function phylandml(
    tree_in::PhyloTree;
    delimiter::String = "_",
    index::Union{Nothing,Int} = nothing,
    regex::Union{Nothing,Regex} = nothing,
    design::Union{Nothing,Matrix{Int}} = nothing,
    start_migrate::Float64 = NaN,
    start_Ne::Float64 = NaN,
    minBL::Float64 = 1e-6,
    likelihood::Symbol = :PL2,
    maxiters::Int = 300
)::PhylandFit
    bt = collapse_singletons(tree_in)
    tts = tips2states(bt.tip_label; delimiter=delimiter, index=index, regex=regex)
    demes = tts.demes
    m = tts.m
    bdt = DatedTree(bt; sampleStates=tts.ssts, demes=demes, minEdgeLength=minBL)
    dm = gen_island_model(demes)

    if design === nothing
        design = reshape(collect(1:(m * m)), m, m)
    else
        design = copy(design)
    end
    for k in 1:m
        design[k, k] = 0
    end

    mig_ids = sort!(unique([design[k, l] for k in 1:m for l in 1:m if k != l && design[k, l] > 0]))
    ne_names = ["Ne$k" for k in 1:m]
    mig_names = ["migrate$id" for id in mig_ids]
    pnames = vcat(ne_names, mig_names)

    Ne0 = isnan(start_Ne) ? max(0.1, bdt.maxHeight / (2.0 * m)) : start_Ne
    mig0 = isnan(start_migrate) ? max(0.01, (bdt.maxHeight / 100.0) / m) : start_migrate
    theta0_log = vcat(fill(log(Ne0), m), fill(log(mig0), length(mig_ids)))

    t1 = bdt.maxSampleTime
    t0 = t1 - bdt.maxHeight - 1.0

    function unpack_theta(log_par::AbstractVector{<:Real})
        d = Dict{String,Float64}()
        for k in 1:m
            d["Ne$k"] = exp(clamp(log_par[k], -15.0, 15.0))
        end
        for k in 1:m, l in 1:m
            if k != l
                id = design[k, l]
                if id <= 0
                    d["mu$(k)$(l)"] = 0.0
                else
                    pos = findfirst(==(id), mig_ids)
                    d["mu$(k)$(l)"] = exp(clamp(log_par[m + pos], -15.0, 10.0))
                end
            end
        end
        return d
    end

    function neg_ll(log_par::Vector{Float64})
        th = unpack_theta(log_par)
        tfgy = dm(th, nothing, t0, t1; res=10)
        ll = colik_pik_fgy(
            bdt, tfgy;
            timeOfOriginBoundaryCondition=true,
            maxHeight=Inf,
            forgiveAgtY=1.0,
            AgtY_penalty=0.0,
            returnTree=false,
            step_size_res=10,
            likelihood=likelihood
        )
        return isfinite(ll) ? -ll : 1e12
    end

    opt = optimize(neg_ll, theta0_log, NelderMead(), Optim.Options(iterations=maxiters))
    best_log = Optim.minimizer(opt)
    best_ll = -Optim.minimum(opt)
    best_theta = unpack_theta(best_log)

    # Ancestral state reconstruction via phydynR::ace
    bdt_ace = ace(bdt, best_theta, dm; t0=t0, res=100, forgiveAgtY=1.0, AgtY_penalty=0.0, step_size_res=10, likelihood=likelihood)
    ace_df = DataFrame(transpose(bdt_ace.acestates), Symbol.(demes))

    # Expected directional lineage transitions across edges from ace posterior states
    exp_jumps = zeros(Float64, m, m)
    for e in 1:size(bdt_ace.tree.edge, 1)
        p = bdt_ace.tree.edge[e, 1]
        c = bdt_ace.tree.edge[e, 2]
        pp = bdt_ace.acestates[:, p]
        pc = bdt_ace.acestates[:, c]
        if all(!isnan, pp) && all(!isnan, pc)
            for k in 1:m, l in 1:m
                k != l && (exp_jumps[k, l] += pp[k] * pc[l])
            end
        end
    end

    nice_coef = Dict{String,Float64}()
    for k in 1:m
        nice_coef[demes[k]] = best_theta["Ne$k"]
    end
    for k in 1:m, l in 1:m
        if k != l
            nice_coef["$(demes[l])->$(demes[k])"] = best_theta["mu$(k)$(l)"]
        end
    end

    return PhylandFit(nice_coef, best_log, pnames, best_ll, demes, design, ace_df, bdt_ace, exp_jumps)
end

"""
    confint_phylandml(fit::PhylandFit, whichparm::Int; guess_se::Float64=1.5, np::Int=13) -> PhylandProfile

Compute profile likelihood and 95% confidence interval for parameter index `whichparm` in `fit::PhylandFit`
(ported from `phyland::confint.phylandml`).
"""
function confint_phylandml(
    fit::PhylandFit,
    whichparm::Int;
    guess_se::Float64 = 1.5,
    np::Int = 13
)::PhylandProfile
    pname = fit.param_names[whichparm]
    center = fit.logcoef[whichparm]
    grid = sort!(unique(vcat(collect(range(center - guess_se, center + guess_se; length=np)), center)))
    dm = gen_island_model(fit.demes)
    m = length(fit.demes)
    mig_ids = sort!(unique([fit.design[k, l] for k in 1:m for l in 1:m if k != l && fit.design[k, l] > 0]))
    t1 = fit.bdt.maxSampleTime
    t0 = t1 - fit.bdt.maxHeight - 1.0

    function eval_at(val::Float64)
        log_par = copy(fit.logcoef)
        log_par[whichparm] = val
        d = Dict{String,Float64}()
        for k in 1:m
            d["Ne$k"] = exp(clamp(log_par[k], -15.0, 15.0))
        end
        for k in 1:m, l in 1:m
            if k != l
                id = fit.design[k, l]
                d["mu$(k)$(l)"] = id <= 0 ? 0.0 : exp(clamp(log_par[m + findfirst(==(id), mig_ids)], -15.0, 10.0))
            end
        end
        tfgy = dm(d, nothing, t0, t1; res=10)
        return colik_pik_fgy(fit.bdt, tfgy; forgiveAgtY=1.0, AgtY_penalty=0.0, step_size_res=10)
    end

    ll = [eval_at(g) for g in grid]
    thresh = maximum(ll) - 1.92
    above = findall(x -> x >= thresh, ll)
    ci_lb = isempty(above) ? NaN : exp(grid[first(above)])
    ci_ub = isempty(above) ? NaN : exp(grid[last(above)])
    sumdat = DataFrame(Parameter_value = exp.(grid), Profile_log_likelihood = ll)
    return PhylandProfile(pname, (ci_lb, ci_ub), grid, ll, sumdat)
end

"""
    boot_phylandml(trees::Vector{PhyloTree}; kwargs...)

Bootstrap `phylandml` across a collection of trees (ported from `phyland::boot.phylandml`).
"""
function boot_phylandml(trees::Vector{PhyloTree}; kwargs...)
    fits = [phylandml(tr; kwargs...) for tr in trees]
    keys_sorted = sort!(collect(keys(fits[1].coef)))
    rows = NamedTuple[]
    for k in keys_sorted
        vals = [f.coef[k] for f in fits]
        push!(rows, (parameter = k, mean = mean(vals), q025 = quantile(vals, 0.025), q975 = quantile(vals, 0.975)))
    end
    return (ci = DataFrame(rows), fits = fits)
end

"""Alias for [`phylandml`](@ref), maximum-likelihood island-model inference."""
const island_coalescent_ml = phylandml
"""Alias for [`confint_phylandml`](@ref), an island-model likelihood interval."""
const confint_island_coalescent = confint_phylandml
"""Alias for [`boot_phylandml`](@ref), island-model inference across replicate trees."""
const boot_island_coalescent = boot_phylandml
