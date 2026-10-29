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

"""
    zoonosis_demographic_process(base_params::ZoonosisParams = model1_params()) -> Function

Return a `phydynR` demographic process function `dm(theta, x0, t0, t1; res=1000)` for the
two-host zoonosis model.
"""
function zoonosis_demographic_process(base_params::ZoonosisParams = model1_params())
    base_mod = ModelFGY(base_params)
    raw_dm = build_demographic_process(base_mod)
    function dm_zoonosis(
        theta::AbstractDict = Dict{String,Float64}(),
        x0::Union{Nothing,AbstractDict} = nothing,
        t0::Real = base_mod.t0,
        t1::Real = base_mod.tfin;
        res::Int = 1000,
        solver = Tsit5()
    )::TFGY
        theta_ode = Dict{String,Float64}(String(k) => Float64(v) for (k, v) in theta)
        psi_val = get(theta_ode, "psi", base_params.psi)
        theta_ode["psi"] = 0.0
        tfgy = raw_dm(theta_ode, x0, t0, t1; res=res, solver=solver)
        if psi_val > 0.0
            mod_tmp = deepcopy(base_mod)
            mod_tmp.t0 = Float64(t0)
            mod_tmp.tfin = Float64(t1)
            for (k, v) in theta_ode
                mod_tmp.parameters[k] = v
            end
            if x0 !== nothing
                for (k, v) in x0
                    mod_tmp.initial[String(k)] = Float64(v)
                end
            end
            sol = solveodes(mod_tmp; solver=solver, res=res)
            beta_val = get(mod_tmp.parameters, "beta", base_params.beta)
            nt = length(tfgy.times)
            for idx in 1:nt
                u = sol.u[nt - idx + 1]
                X1_val = max(0.0, u[3])
                Y2_val = max(0.0, u[2])
                tfgy.Fs[idx][2, 1] = psi_val * beta_val * X1_val * Y2_val
            end
        end
        return tfgy
    end
    return dm_zoonosis
end
