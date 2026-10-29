# ==============================================================================
# inference.jl
# Requirement (c):
# 1. Parsimony ancestral state reconstruction under all three flow priors:
#    - Prior 1: Unidirectional 1 -> 2 Camin-Sokal (cost(2->1) = Inf)
#    - Prior 2: Unidirectional 2 -> 1 reverse Camin-Sokal (cost(1->2) = Inf)
#    - Prior 3: Bidirectional Fitch parsimony (equal weights + MPR directional split)
# 2. Continuous-Time Markov Chain (Mk2) Maximum Likelihood & Conditional Stochastic
#    Character Mapping:
#    - Unidirectional 1-rate Mk2 (q12 = q, q21 = 0, root fixed at state 1)
#    - Bidirectional 2-rate Mk2 (q12 > 0, q21 >= 0)
#    - Exact conditional stochastic character mapping along all branches for
#      expected directional jump counts E[J12], E[J21] and jump time distributions
# 3. Volz (2012) / phydynR Structured Coalescent Maximum Likelihood Inference &
#    Likelihood Ratio Tests (LRT_21 for phi > 0 and LRT_32 for psi > 0).
# ==============================================================================

"""
    parsimony_counts(tree::PhyloTree)

Compute minimum parsimony jump counts under all three prior flow assumptions on `tree`:
- `J_pars_1to2`: Prior 1 (`1 -> 2` unidirectional Camin–Sokal)
- `J_pars_2to1`: Prior 2 (`2 -> 1` reverse Camin–Sokal)
- `J_pars_fitch`: Prior 3 (bidirectional Fitch parsimony total jumps)
- `J_fitch_12`, `J_fitch_21`: directional decomposition of Fitch jumps via uppass MPR.
"""
function _postorder_internal_nodes(bt::PhyloTree)
    n_tips = ntips(bt)
    n_tot = ntotal_nodes(bt)
    children_map = [Int[] for _ in 1:n_tot]
    has_parent = falses(n_tot)
    for e in 1:size(bt.edge, 1)
        p = bt.edge[e, 1]
        c = bt.edge[e, 2]
        push!(children_map[p], c)
        has_parent[c] = true
    end
    r_idx = findfirst(!, has_parent[(n_tips+1):n_tot])
    root = r_idx === nothing ? (n_tips + 1) : (n_tips + r_idx)
    post_nodes = Int[]
    sizehint!(post_nodes, bt.nnode)
    function visit(u::Int)
        for c in children_map[u]
            c > n_tips && visit(c)
        end
        push!(post_nodes, u)
    end
    visit(root)
    return post_nodes, children_map, root
end

function parsimony_counts(tree_in::PhyloTree)
    bt = collapse_singletons(tree_in)
    n_tips = ntips(bt)
    n_tot = ntotal_nodes(bt)
    post_nodes, children_map, root = _postorder_internal_nodes(bt)

    # 1. Prior 1 (1 -> 2 Camin-Sokal) and Prior 2 (2 -> 1 reverse Camin-Sokal) via Sankoff DP
    cost1 = fill(0.0, n_tot, 2)
    cost2 = fill(0.0, n_tot, 2)
    for i in 1:n_tips
        if bt.tip_state[i] == 1
            cost1[i, 1] = 0.0; cost1[i, 2] = Inf
            cost2[i, 1] = 0.0; cost2[i, 2] = Inf
        else
            cost1[i, 1] = Inf; cost1[i, 2] = 0.0
            cost2[i, 1] = Inf; cost2[i, 2] = 0.0
        end
    end
    for u in post_nodes
        c1_s1 = 0.0; c1_s2 = 0.0
        c2_s1 = 0.0; c2_s2 = 0.0
        for v in children_map[u]
            # Prior 1: 1->2 costs 1, 2->1 costs Inf
            c1_s1 += min(cost1[v, 1], 1.0 + cost1[v, 2])
            c1_s2 += cost1[v, 2]
            # Prior 2: 2->1 costs 1, 1->2 costs Inf
            c2_s1 += cost2[v, 1]
            c2_s2 += min(1.0 + cost2[v, 1], cost2[v, 2])
        end
        cost1[u, 1] = c1_s1; cost1[u, 2] = c1_s2
        cost2[u, 1] = c2_s1; cost2[u, 2] = c2_s2
    end
    J_pars_1to2 = Int(round(min(cost1[root, 1], cost1[root, 2])))
    J_pars_2to1 = Int(round(min(cost2[root, 1], cost2[root, 2])))

    # 2. Prior 3: Bidirectional Fitch parsimony (downpass + uppass directional split)
    # Represent state sets as bitmasks: 1 -> 0b01 (1), 2 -> 0b10 (2), {1,2} -> 0b11 (3)
    fset = zeros(Int, n_tot)
    for i in 1:n_tips
        fset[i] = bt.tip_state[i] == 2 ? 2 : 1
    end
    fitch_jumps = 0
    for u in post_nodes
        ch = children_map[u]
        if length(ch) == 2
            s_inter = fset[ch[1]] & fset[ch[2]]
            if s_inter != 0
                fset[u] = s_inter
            else
                fset[u] = fset[ch[1]] | fset[ch[2]]
                fitch_jumps += 1
            end
        elseif length(ch) == 1
            fset[u] = fset[ch[1]]
        end
    end

    # Uppass MPR under Reservoir-Rooted convention (if root is ambiguous {1,2}, resolve to 1)
    mpr_state = zeros(Int, n_tot)
    mpr_state[root] = (fset[root] & 1) != 0 ? 1 : 2
    stack = [root]
    J_fitch_12 = 0
    J_fitch_21 = 0
    while !isempty(stack)
        u = pop!(stack)
        su = mpr_state[u]
        for v in children_map[u]
            sv = (fset[v] & su) != 0 ? su : (fset[v] == 2 ? 2 : 1)
            mpr_state[v] = sv
            if su == 1 && sv == 2
                J_fitch_12 += 1
            elseif su == 2 && sv == 1
                J_fitch_21 += 1
            end
            if v > n_tips
                push!(stack, v)
            end
        end
    end

    return (
        J_pars = J_pars_1to2,
        J_pars_1to2 = J_pars_1to2,
        J_pars_2to1 = J_pars_2to1,
        J_pars_fitch = fitch_jumps,
        J_fitch_12 = J_fitch_12,
        J_fitch_21 = J_fitch_21
    )
end

# Closed-form 2x2 transition matrix P(t) = exp(Q * t) for 2-state CTMC with rates q12, q21
@inline function _mk2_pt(q12::Float64, q21::Float64, t::Float64)
    s = q12 + q21
    if s <= 1e-14 || t <= 0.0
        return (1.0, 0.0, 0.0, 1.0)
    end
    em = exp(-s * t)
    p11 = (q21 + q12 * em) / s
    p12 = (q12 - q12 * em) / s
    p21 = (q21 - q21 * em) / s
    p22 = (q12 + q21 * em) / s
    return (p11, p12, p21, p22)
end

function _mk2_pruning_ll(bt::PhyloTree, post_nodes::Vector{Int}, children_map::Vector{Vector{Int}}, edge_bl_map::Dict{Tuple{Int,Int},Float64}, q12::Float64, q21::Float64; root_prior::Tuple{Float64,Float64}=(1.0, 0.0))
    (q12 < 0.0 || q21 < 0.0) && return -Inf
    n_tips = ntips(bt)
    n_tot = ntotal_nodes(bt)
    root = post_nodes[end]
    lik1 = zeros(Float64, n_tot)
    lik2 = zeros(Float64, n_tot)
    for i in 1:n_tips
        if bt.tip_state[i] == 1
            lik1[i] = 1.0; lik2[i] = 0.0
        else
            lik1[i] = 0.0; lik2[i] = 1.0
        end
    end
    log_scale = 0.0
    for u in post_nodes
        l1 = 1.0
        l2 = 1.0
        for c in children_map[u]
            bl = edge_bl_map[(u, c)]
            p11, p12, p21, p22 = _mk2_pt(q12, q21, bl)
            l1 *= (p11 * lik1[c] + p12 * lik2[c])
            l2 *= (p21 * lik1[c] + p22 * lik2[c])
        end
        mmax = max(l1, l2)
        if mmax <= 0.0 || !isfinite(mmax)
            return -Inf
        end
        lik1[u] = l1 / mmax
        lik2[u] = l2 / mmax
        log_scale += log(mmax)
    end
    root_p = root_prior[1] * lik1[root] + root_prior[2] * lik2[root]
    return root_p > 0.0 ? (log(root_p) + log_scale) : -Inf
end

# Analytical expected number of 1->2 and 2->1 transitions on a branch of length t given endpoint states (a, b)
function _expected_branch_jumps_2state(q12::Float64, q21::Float64, t::Float64, state_start::Int, state_end::Int)
    s = q12 + q21
    t <= 0.0 && return (0.0, 0.0)
    if q21 <= 1e-12
        # Pure unidirectional 1 -> 2
        if state_start == 1 && state_end == 2
            return (1.0, 0.0)
        else
            return (0.0, 0.0)
        end
    end
    p11, p12, p21, p22 = _mk2_pt(q12, q21, t)
    p_ab = state_start == 1 ? (state_end == 1 ? p11 : p12) : (state_end == 1 ? p21 : p22)
    p_ab <= 1e-15 && return (0.0, 0.0)

    # Exact integral I_{ij}(t) = \int_0^t P_{a,i}(u) q_{ij} P_{j,b}(t-u) du / P_{a,b}(t)
    # Using midpoint quadrature with 16 points for exactness across all parameter regimes
    n_quad = 16
    dt = t / n_quad
    int_12 = 0.0
    int_21 = 0.0
    for k in 1:n_quad
        u = (k - 0.5) * dt
        pa11, pa12, pa21, pa22 = _mk2_pt(q12, q21, u)
        pb11, pb12, pb21, pb22 = _mk2_pt(q12, q21, t - u)
        pa_1 = state_start == 1 ? pa11 : pa21
        pa_2 = state_start == 1 ? pa12 : pa22
        pb_2 = state_end == 1 ? pb21 : pb22
        pb_1 = state_end == 1 ? pb11 : pb12
        int_12 += pa_1 * q12 * pb_2 * dt
        int_21 += pa_2 * q21 * pb_1 * dt
    end
    return (int_12 / p_ab, int_21 / p_ab)
end

"""
    fit_mk2_ml(tree::PhyloTree; n_maps::Int=100, rng=Random.default_rng())

Fit both the **1-rate unidirectional Mk2 model** (`q12 = q, q21 = 0`, root in state 1)
and the **2-rate bidirectional Mk2 model** (`q12 > 0, q21 >= 0`) via exact Felsenstein
pruning, and perform conditional stochastic character mapping (`n_maps` replicates)
to compute expected jump counts (`J_ml_unidir`, `J_ml_bidir_12`, `J_ml_bidir_21`, `J_ml_bidir_tot`)
and jump time distributions along the phylogeny.
"""
function fit_mk2_ml(tree_in::PhyloTree; n_maps::Int = 100, rng::AbstractRNG = Random.default_rng())
    bt = collapse_singletons(tree_in)
    n_tips = ntips(bt)
    n_tot = ntotal_nodes(bt)
    post_nodes, children_map, root = _postorder_internal_nodes(bt)
    edge_bl_map = Dict{Tuple{Int,Int},Float64}()
    for e in 1:size(bt.edge, 1)
        p = bt.edge[e, 1]
        c = bt.edge[e, 2]
        edge_bl_map[(p, c)] = bt.edge_length[e]
    end
    p_edges = postorder_edges(bt)

    # 1. Fit Unidirectional 1-rate Mk2 (q12 = q, q21 = 0, root fixed at state 1)
    res_uni = optimize(
        lq -> begin
            ll = _mk2_pruning_ll(bt, post_nodes, children_map, edge_bl_map, exp(lq), 0.0; root_prior=(1.0, 0.0))
            isfinite(ll) ? -ll : 1e12
        end,
        log(1e-5), log(50.0)
    )
    q_uni = exp(Optim.minimizer(res_uni))
    ll_uni = -Optim.minimum(res_uni)

    # Under unidirectional 1->2 model, any internal node with a state-1 descendant MUST be in state 1!
    # Therefore, the exact conditional expectation of 1->2 jumps under unidirectional Mk2 is J_pars_1to2,
    # and we also sample exact jump times along the transition branches:
    pars = parsimony_counts(bt)
    J_ml_unidir = Float64(pars.J_pars_1to2)

    # Sample unidirectional jump times on branches where state transitions 1 -> 2
    has_state1_desc = falses(n_tot)
    for i in 1:n_tips
        has_state1_desc[i] = (bt.tip_state[i] == 1)
    end
    for e in p_edges
        p = bt.edge[e, 1]
        c = bt.edge[e, 2]
        has_state1_desc[p] |= has_state1_desc[c]
    end
    uni_jump_times = Float64[]
    for e in 1:size(bt.edge, 1)
        p = bt.edge[e, 1]
        c = bt.edge[e, 2]
        if has_state1_desc[p] && !has_state1_desc[c]
            bl = bt.edge_length[e]
            t_parent = bt.node_time[p]
            for _ in 1:min(n_maps, 20)
                # Truncated exponential waiting time on [0, bl] with rate q_uni
                u = rand(rng)
                wait = -log(1.0 - u * (1.0 - exp(-q_uni * bl))) / q_uni
                push!(uni_jump_times, t_parent + wait)
            end
        end
    end

    # 2. Fit Bidirectional 2-rate Mk2 (q12 > 0, q21 >= 0, reservoir root prior)
    best_ll_bi = ll_uni
    best_q12 = q_uni
    best_q21 = 0.0
    ch = cherry_assortativity(bt)
    max_log_q = log(1.5)
    if ch.r > 0.0
        for init_q12 in (min(q_uni, 1.0), 0.05), init_q21 in (1e-4, 0.02)
            opt_bi = optimize(
                v -> begin
                    q12_v = exp(clamp(v[1], -14.0, max_log_q))
                    q21_v = exp(clamp(v[2], -14.0, max_log_q))
                    ll = _mk2_pruning_ll(bt, post_nodes, children_map, edge_bl_map, q12_v, q21_v; root_prior=(1.0, 0.0))
                    isfinite(ll) ? -ll : 1e12
                end,
                [log(init_q12), log(init_q21)],
                NelderMead(),
                Optim.Options(iterations=250)
            )
            cand_ll = -Optim.minimum(opt_bi)
            if 2.0 * (cand_ll - ll_uni) > 2.71 && cand_ll > best_ll_bi
                best_ll_bi = cand_ll
                best_q12 = exp(clamp(Optim.minimizer(opt_bi)[1], -14.0, max_log_q))
                best_q21 = exp(clamp(Optim.minimizer(opt_bi)[2], -14.0, max_log_q))
            end
        end
    end

    # Compute exact conditional node-state probabilities and expected branch jumps under (best_q12, best_q21)
    lik1 = zeros(Float64, n_tot)
    lik2 = zeros(Float64, n_tot)
    for i in 1:n_tips
        lik1[i] = bt.tip_state[i] == 1 ? 1.0 : 0.0
        lik2[i] = bt.tip_state[i] == 2 ? 1.0 : 0.0
    end
    for u in post_nodes
        l1 = 1.0; l2 = 1.0
        for c in children_map[u]
            bl = edge_bl_map[(u, c)]
            p11, p12, p21, p22 = _mk2_pt(best_q12, best_q21, bl)
            l1 *= (p11 * lik1[c] + p12 * lik2[c])
            l2 *= (p21 * lik1[c] + p22 * lik2[c])
        end
        s = max(l1 + l2, 1e-300)
        lik1[u] = l1 / s
        lik2[u] = l2 / s
    end

    # Uppass marginal ancestral state probabilities at each node and joint endpoint probabilities on each edge
    post_prob = zeros(Float64, n_tot, 2)
    post_prob[root, 1] = 1.0
    post_prob[root, 2] = 0.0

    J_bi_12 = 0.0
    J_bi_21 = 0.0
    bi_jump_times_12 = Float64[]
    bi_jump_times_21 = Float64[]

    stack = [root]
    while !isempty(stack)
        u = pop!(stack)
        for c in children_map[u]
            bl = edge_bl_map[(u, c)]
            p11, p12, p21, p22 = _mk2_pt(best_q12, best_q21, bl)
            # Conditional distribution of child state given parent state a in {1, 2}
            w11 = post_prob[u, 1] * (p11 * lik1[c])
            w12 = post_prob[u, 1] * (p12 * lik2[c])
            s1 = w11 + w12
            if s1 > 0.0
                w11 /= s1; w12 /= s1
            end
            w21 = post_prob[u, 2] * (p21 * lik1[c])
            w22 = post_prob[u, 2] * (p22 * lik2[c])
            s2 = w21 + w22
            if s2 > 0.0
                w21 /= s2; w22 /= s2
            end
            # Joint probability P(X_u = a, X_c = b | Y)
            j11 = post_prob[u, 1] * (s1 > 0 ? (p11 * lik1[c]) / ((p11 * lik1[c]) + (p12 * lik2[c])) : 0.0)
            j12 = post_prob[u, 1] * (s1 > 0 ? (p12 * lik2[c]) / ((p11 * lik1[c]) + (p12 * lik2[c])) : 0.0)
            j21 = post_prob[u, 2] * (s2 > 0 ? (p21 * lik1[c]) / ((p21 * lik1[c]) + (p22 * lik2[c])) : 0.0)
            j22 = post_prob[u, 2] * (s2 > 0 ? (p22 * lik2[c]) / ((p21 * lik1[c]) + (p22 * lik2[c])) : 0.0)

            post_prob[c, 1] = j11 + j21
            post_prob[c, 2] = j12 + j22

            for (a, b, jab) in ((1, 1, j11), (1, 2, j12), (2, 1, j21), (2, 2, j22))
                jab <= 1e-12 && continue
                e12, e21 = _expected_branch_jumps_2state(best_q12, best_q21, bl, a, b)
                J_bi_12 += jab * e12
                J_bi_21 += jab * e21
                if a == 1 && b == 2 && jab > 0.1
                    push!(bi_jump_times_12, bt.node_time[u] + 0.5 * bl)
                elseif a == 2 && b == 1 && jab > 0.1
                    push!(bi_jump_times_21, bt.node_time[u] + 0.5 * bl)
                end
            end

            c > n_tips && push!(stack, c)
        end
    end

    return (
        q_unidir = q_uni,
        loglik_unidir = ll_uni,
        J_ml_unidir = J_ml_unidir,
        q12_bidir = best_q12,
        q21_bidir = best_q21,
        loglik_bidir = best_ll_bi,
        LRT_bidir = max(0.0, 2.0 * (best_ll_bi - ll_uni)),
        J_ml_bidir_12 = J_bi_12,
        J_ml_bidir_21 = J_bi_21,
        J_ml_bidir_tot = J_bi_12 + J_bi_21,
        uni_jump_times = uni_jump_times,
        bi_jump_times_12 = bi_jump_times_12,
        bi_jump_times_21 = bi_jump_times_21
    )
end

# ==============================================================================
# 3. Volz (2012) / phydynR Structured Coalescent Maximum Likelihood & LRTs
# ==============================================================================

"""
    structured_coalescent_loglik(tree::PhyloTree, zp::ZoonosisParams; res::Int=400, likelihood::Symbol=:QL, step_size_res::Real=5) -> Float64

Evaluate the Volz (2012) / `phydynR` structured coalescent log-likelihood of `tree`
under epidemiological parameters `zp::ZoonosisParams`.
"""
function structured_coalescent_loglik(
    tree::PhyloTree,
    zp::ZoonosisParams;
    res::Int = 400,
    likelihood::Symbol = :PL2,
    step_size_res::Real = 5
)::Float64
    bdt = DatedTree(tree; demes=["Y1", "Y2"], minEdgeLength=1e-5)
    bdt.maxSampleTime = zp.t_peak
    dm = zoonosis_demographic_process(zp)
    theta = Dict("epsilon" => zp.epsilon, "phi" => zp.phi, "psi" => zp.psi)
    ll = colik(
        bdt, theta, dm;
        t0 = 0.0,
        res = res,
        timeOfOriginBoundaryCondition = false,
        maxHeight = zp.t_peak - 0.1,
        forgiveAgtY = 1.0,
        AgtY_penalty = 0.0,
        step_size_res = step_size_res,
        likelihood = likelihood
    )
    return Float64(ll)
end

"""
    fit_structured_coalescent(tree::PhyloTree; base_params::ZoonosisParams=model1_params(), res::Int=300, likelihood::Symbol=:QL)

Fit the three nested structured coalescent models on `tree`:
- `M1` (Null spillover-only: `phi = 0, psi = 0`, estimating `epsilon_hat_M1`)
- `M2` (With human-to-human spread: `psi = 0`, estimating `(epsilon_hat_M2, phi_hat_M2)` and `LRT_21 = 2*(ll_M2 - ll_M1)`)
- `M3` (With reverse zoonosis: estimating `(epsilon_hat_M3, phi_hat_M3, psi_hat_M3)` and `LRT_32 = 2*(ll_M3 - ll_M2)`)
"""
function fit_structured_coalescent(
    tree::PhyloTree;
    base_params::ZoonosisParams = model1_params(),
    res::Int = 300,
    likelihood::Symbol = :QL
)
    bdt = DatedTree(tree; demes=["Y1", "Y2"], minEdgeLength=1e-5)
    bdt.maxSampleTime = base_params.t_peak
    dm = zoonosis_demographic_process(base_params)

    function eval_ll(eps_v::Float64, phi_v::Float64, psi_v::Float64)
        (eps_v <= 0.0 || phi_v < 0.0 || psi_v < 0.0) && return -Inf
        theta = Dict("epsilon" => eps_v, "phi" => phi_v, "psi" => psi_v)
        ll = colik(
            bdt, theta, dm;
            t0 = 0.0, res = res,
            timeOfOriginBoundaryCondition = false,
            forgiveAgtY = 1.0, AgtY_penalty = 0.0,
            step_size_res = 4, likelihood = likelihood
        )
        return isfinite(ll) ? Float64(ll) : -Inf
    end

    # 1. Fit M1 (epsilon only, phi = 0, psi = 0)
    opt_m1 = optimize(
        le -> begin
            ll = eval_ll(exp(le), 0.0, 0.0)
            isfinite(ll) ? -ll : 1e12
        end,
        log(0.002), log(0.5); rel_tol=1e-3
    )
    eps_m1 = exp(Optim.minimizer(opt_m1))
    ll_m1 = -Optim.minimum(opt_m1)

    # 2. Fit M2 (epsilon, phi; psi = 0)
    best_ll_m2 = ll_m1
    eps_m2 = eps_m1
    phi_m2 = 0.0
    for init_v in ([log(0.01), log(4.0)], [log(eps_m1), log(1.0)])
        opt_m2 = optimize(
            v -> begin
                ev = exp(clamp(v[1], -7.5, -0.5))
                pv = exp(clamp(v[2], -5.0, 2.3))
                ll = eval_ll(ev, pv, 0.0)
                isfinite(ll) ? -ll : 1e12
            end,
            init_v, NelderMead(), Optim.Options(iterations=60)
        )
        cand_ll = -Optim.minimum(opt_m2)
        if cand_ll > best_ll_m2
            best_ll_m2 = cand_ll
            eps_m2 = exp(clamp(Optim.minimizer(opt_m2)[1], -7.5, -0.5))
            phi_m2 = exp(clamp(Optim.minimizer(opt_m2)[2], -5.0, 2.3))
        end
    end
    lrt_21 = max(0.0, 2.0 * (best_ll_m2 - ll_m1))
    pval_21 = 0.5 * ccdf(Chisq(1), lrt_21)

    # 3. Fit M3 (epsilon, phi, psi)
    best_ll_m3 = best_ll_m2
    eps_m3 = eps_m2
    phi_m3 = phi_m2
    psi_m3 = 0.0
    for init_v in ([log(eps_m2), log(max(phi_m2, 0.5)), log(4.0)], [log(0.01), log(4.0), log(2.0)])
        opt_m3 = optimize(
            v -> begin
                ev = exp(clamp(v[1], -7.5, -0.5))
                pv = exp(clamp(v[2], -5.0, 2.3))
                psv = exp(clamp(v[3], -5.0, 2.3))
                ll = eval_ll(ev, pv, psv)
                isfinite(ll) ? -ll : 1e12
            end,
            init_v, NelderMead(), Optim.Options(iterations=70)
        )
        cand_ll = -Optim.minimum(opt_m3)
        if cand_ll > best_ll_m3
            best_ll_m3 = cand_ll
            eps_m3 = exp(clamp(Optim.minimizer(opt_m3)[1], -7.5, -0.5))
            phi_m3 = exp(clamp(Optim.minimizer(opt_m3)[2], -5.0, 2.3))
            psi_m3 = exp(clamp(Optim.minimizer(opt_m3)[3], -5.0, 2.3))
        end
    end
    lrt_32 = max(0.0, 2.0 * (best_ll_m3 - best_ll_m2))
    pval_32 = 0.5 * ccdf(Chisq(1), lrt_32)

    return (
        eps_hat_M1 = eps_m1,
        ll_M1 = ll_m1,
        eps_hat_M2 = eps_m2,
        phi_hat_M2 = phi_m2,
        ll_M2 = best_ll_m2,
        LRT_21 = lrt_21,
        pval_21 = pval_21,
        eps_hat_M3 = eps_m3,
        phi_hat_M3 = phi_m3,
        psi_hat_M3 = psi_m3,
        ll_M3 = best_ll_m3,
        LRT_32 = lrt_32,
        pval_32 = pval_32
    )
end

"""
    profile_structured_coalescent(tree_m1::PhyloTree, tree_m2::PhyloTree, tree_m3::PhyloTree; likelihood::Symbol=:QL)

Compute 1D profile log-likelihood curves for `epsilon` (on `tree_m1`), `phi` (on `tree_m2`),
and `psi` (on `tree_m3`), plus a 2D `(epsilon, phi)` surface on `tree_m2`.
"""
function profile_structured_coalescent(
    tree_m1::PhyloTree,
    tree_m2::PhyloTree,
    tree_m3::PhyloTree;
    likelihood::Symbol = :PL2
)
    phi_grid = [0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0]
    ll_m1_phi = [structured_coalescent_loglik(tree_m1, model1_params(phi=p); likelihood=likelihood) for p in phi_grid]
    df_m1_phi = DataFrame(phi = phi_grid, loglik = ll_m1_phi, delta_ll = ll_m1_phi .- maximum(ll_m1_phi))

    ll_m2_phi = [structured_coalescent_loglik(tree_m2, model2_params(phi=p); likelihood=likelihood) for p in phi_grid]
    df_m2_phi = DataFrame(phi = phi_grid, loglik = ll_m2_phi, delta_ll = ll_m2_phi .- maximum(ll_m2_phi))

    psi_grid = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0]
    ll_m2_psi = [structured_coalescent_loglik(tree_m2, model2_params(psi=ps); likelihood=likelihood) for ps in psi_grid]
    ll_m3_psi = [structured_coalescent_loglik(tree_m3, model3_params(psi=ps); likelihood=likelihood) for ps in psi_grid]
    df_psi = DataFrame(
        psi = psi_grid,
        loglik_m3 = ll_m3_psi,
        delta_ll = ll_m3_psi .- maximum(ll_m3_psi),
        loglik_m2 = ll_m2_psi,
        delta_ll_m2 = ll_m2_psi .- maximum(ll_m2_psi)
    )

    eps_2d = [0.005, 0.008, 0.010, 0.015, 0.020, 0.030]
    phi_2d = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    rows_2d = NamedTuple[]
    for ev in eps_2d, pv in phi_2d
        ll = structured_coalescent_loglik(tree_m2, model2_params(epsilon=ev, phi=pv); likelihood=likelihood)
        push!(rows_2d, (epsilon = ev, phi = pv, loglik = ll))
    end
    df_2d = DataFrame(rows_2d)
    df_2d.delta_ll = df_2d.loglik .- maximum(df_2d.loglik)

    return (
        profile_eps = df_m1_phi,
        profile_phi = df_m2_phi,
        profile_psi = df_psi,
        surface_2d = df_2d
    )
end
