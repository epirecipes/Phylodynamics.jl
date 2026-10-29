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
