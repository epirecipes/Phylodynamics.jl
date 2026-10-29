"""
    build_zoonosis_reaction_system(; name::Symbol=:zoonosis)

Construct the two-host SIRS zoonosis model (covering Models 1, 2, and 3) using
`ModelingToolkit.jl` / `Catalyst.jl` reaction system syntax (`@independent_variables`, `@parameters`,
`@species`, `Reaction`, `ReactionSystem`), matching Simon Frost's `rn_mtk.qmd` style.

Returns `(rs = rs, species = (X1=X1, Y1=Y1, X2=X2, Y2=Y2), params = (...), t = t)`.
"""
function build_zoonosis_reaction_system(; name::Symbol = :zoonosis)
    t = ModelingToolkit.t_nounits
    @parameters nu1 mu1 gamma1 m1 beta nu2 mu2 gamma2 m2 epsilon phi psi
    @species X1(t) Y1(t) X2(t) Y2(t)

    rxs = [
        # 1. Reservoir demographic turnover
        Catalyst.Reaction(nu1, nothing, [X1], nothing, [1]),
        Catalyst.Reaction(mu1, [X1], nothing, [1], nothing),
        # 2. Reservoir-to-reservoir transmission (1 -> 1): X1 + Y1 --> 2Y1
        Catalyst.Reaction(beta, [X1, Y1], [Y1], [1, 1], [2]),
        # 3. Reservoir recovery & mortality: Y1 --> 0
        Catalyst.Reaction(gamma1, [Y1], nothing, [1], nothing),
        Catalyst.Reaction(mu1 + m1, [Y1], nothing, [1], nothing),
        # 4. Human demographic turnover
        Catalyst.Reaction(nu2, nothing, [X2], nothing, [1]),
        Catalyst.Reaction(mu2, [X2], nothing, [1], nothing),
        # 5. Zoonotic spillover (1 -> 2): X2 + Y1 --> Y2 + Y1
        Catalyst.Reaction(epsilon * beta, [X2, Y1], [Y2, Y1], [1, 1], [1, 1]),
        # 6. Human-to-human transmission (2 -> 2): X2 + Y2 --> 2Y2
        Catalyst.Reaction(phi * beta, [X2, Y2], [Y2], [1, 1], [2]),
        # 7. Reverse zoonosis / zooanthroponosis (2 -> 1): X1 + Y2 --> Y1 + Y2
        Catalyst.Reaction(psi * beta, [X1, Y2], [Y1, Y2], [1, 1], [1, 1]),
        # 8. Human recovery & mortality: Y2 --> 0
        Catalyst.Reaction(gamma2, [Y2], nothing, [1], nothing),
        Catalyst.Reaction(mu2 + m2, [Y2], nothing, [1], nothing)
    ]

    rs = Catalyst.ReactionSystem(
        rxs, t, [X1, Y1, X2, Y2],
        [nu1, mu1, gamma1, m1, beta, nu2, mu2, gamma2, m2, epsilon, phi, psi];
        name = name,
        combinatoric_ratelaws = false
    )
    rs = complete(rs)
    return (
        rs = rs,
        species = (X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2),
        infected = [Y1, Y2],
        params = (nu1 = nu1, mu1 = mu1, gamma1 = gamma1, m1 = m1, beta = beta,
                  nu2 = nu2, mu2 = mu2, gamma2 = gamma2, m2 = m2,
                  epsilon = epsilon, phi = phi, psi = psi),
        t = t
    )
end

"""
    build_zoonosis_catalyst_network()

Construct the two-host SIRS zoonosis model using `Catalyst.@reaction_network` DSL.
"""
function build_zoonosis_catalyst_network()
    rn = @reaction_network zoonosis_catalyst begin
        @parameters nu1 mu1 gamma1 m1 beta nu2 mu2 gamma2 m2 epsilon phi psi
        nu1, 0 --> X1
        mu1, X1 --> 0
        beta, X1 + Y1 --> 2Y1
        gamma1, Y1 --> 0
        mu1 + m1, Y1 --> 0
        nu2, 0 --> X2
        mu2, X2 --> 0
        epsilon * beta, X2 + Y1 --> Y2 + Y1
        phi * beta, X2 + Y2 --> 2Y2
        psi * beta, X1 + Y2 --> Y1 + Y2
        gamma2, Y2 --> 0
        mu2 + m2, Y2 --> 0
    end
    return rn
end

function _mtk_u0_and_p(sys_bundle, zp::ZoonosisParams; int_u0::Bool=false)
    sp = sys_bundle.species
    pr = sys_bundle.params
    u0 = int_u0 ? [
        sp.X1 => round(Int, zp.X1_0),
        sp.Y1 => round(Int, zp.Y1_0),
        sp.X2 => round(Int, zp.X2_0),
        sp.Y2 => round(Int, zp.Y2_0)
    ] : [
        sp.X1 => zp.X1_0,
        sp.Y1 => zp.Y1_0,
        sp.X2 => zp.X2_0,
        sp.Y2 => zp.Y2_0
    ]
    p = [
        pr.nu1 => zp.nu1, pr.mu1 => zp.mu1, pr.gamma1 => zp.gamma1, pr.m1 => zp.m1,
        pr.beta => zp.beta, pr.nu2 => zp.nu2, pr.mu2 => zp.mu2, pr.gamma2 => zp.gamma2,
        pr.m2 => zp.m2, pr.epsilon => zp.epsilon, pr.phi => zp.phi, pr.psi => zp.psi
    ]
    return u0, p
end

"""
    simulate_ode_trajectory(zp::ZoonosisParams; tspan=(0.0, zp.t_peak), solver=Tsit5(), saveat=nothing)

Convert the `ModelingToolkit.ReactionSystem` to an `ODESystem` (`Catalyst.ode_model(rs)`)
and solve the forward ODEs using `OrdinaryDiffEqTsit5` (`Tsit5()`).
"""
function simulate_ode_trajectory(
    zp::ZoonosisParams;
    tspan::Tuple{Float64,Float64} = (0.0, zp.t_peak),
    solver = Tsit5(),
    saveat = nothing,
    reltol::Float64 = 1e-9,
    abstol::Float64 = 1e-11
)
    bundle = build_zoonosis_reaction_system()
    odesys = complete(Catalyst.ode_model(bundle.rs))
    u0, p = _mtk_u0_and_p(bundle, zp; int_u0=false)
    oprob = ODEProblem(bundle.rs, u0, tspan, p)
    if saveat === nothing
        sol = solve(oprob, solver; reltol=reltol, abstol=abstol)
    else
        sol = solve(oprob, solver; saveat=saveat, reltol=reltol, abstol=abstol)
    end
    return (sol = sol, bundle = bundle, odesys = odesys)
end

"""
    simulate_stochastic_trajectory(zp::ZoonosisParams; tspan=(0.0, zp.t_peak), rng=Random.default_rng())

Convert the `ModelingToolkit.ReactionSystem` to a `JumpProblem(rs, u0, tspan, p)`
and solve the exact Gillespie stochastic trajectory using `JumpProcesses.jl` (`SSAStepper()`).
"""
function simulate_stochastic_trajectory(
    zp::ZoonosisParams;
    tspan::Tuple{Float64,Float64} = (0.0, zp.t_peak),
    rng::AbstractRNG = Random.default_rng()
)
    bundle = build_zoonosis_reaction_system()
    u0i, p = _mtk_u0_and_p(bundle, zp; int_u0=true)
    jprob = JumpProblem(bundle.rs, u0i, tspan, p; rng=rng)
    jsol = solve(jprob, SSAStepper())
    return (sol = jsol, bundle = bundle, jprob = jprob)
end

"""
    simulate_coalescent_tree(zp::ZoonosisParams; n1::Int=100, n2::Int=100, T_sample::Float64=zp.t_peak, rng=Random.default_rng()) -> PhyloTree

Simulate a structured coalescent tree from `zp::ZoonosisParams` using `Coalescent.jl`'s
`ModelFGY` + `SampleConfiguration` + `SimTree` algorithm and convert to `PhyloTree`.
"""
function simulate_coalescent_tree(
    zp::ZoonosisParams;
    n1::Int = 100,
    n2::Int = 100,
    T_sample::Float64 = zp.t_peak,
    rng::AbstractRNG = Random.default_rng()
)::PhyloTree
    bundle = build_zoonosis_reaction_system()
    u0, p = _mtk_u0_and_p(bundle, zp; int_u0=false)
    sp = bundle.species
    samples = [sp.Y1 => n1, sp.Y2 => n2]
    mod = ModelFGY(bundle.rs, bundle.infected, u0, (0.0, T_sample), p)
    scfg = SampleConfiguration(samples, T_sample)
    st = SimTree(mod, scfg; computedescendants=true, rng=rng)
    return simtree_to_phylotree(st; T_sample=T_sample)
end

"""
    simulate_remaster_tree(zp::ZoonosisParams; n1=100, n2=100, T_sample=zp.t_peak, kwargs...)

Build the two-host study reaction system and delegate to the generic
[`simulate_remaster_tree`](@ref) forward-backward simulator.
"""
function simulate_remaster_tree(zp::ZoonosisParams;
    n1::Int=100, n2::Int=100, T_sample::Float64=zp.t_peak, kwargs...)
    bundle = build_zoonosis_reaction_system()
    u0, p = _mtk_u0_and_p(bundle, zp; int_u0=true)
    return simulate_remaster_tree(bundle.rs, bundle.infected, u0,
        (0.0, T_sample), p, [bundle.species.Y1 => n1, bundle.species.Y2 => n2];
        kwargs...)
end
