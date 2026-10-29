# ==============================================================================
# coalescent_odes.jl
# Analytical Frost/Volz forward epidemiological ODEs and 37-state backward
# structured coalescent ODEs across all three epidemiological models (Models 1, 2, 3)
# and all three parsimony flow priors (Prior 1: 1->2 Camin-Sokal, Prior 2: 2->1
# reverse Camin-Sokal, Prior 3: 6-state bidirectional Fitch), plus cherry
# assortativity r and normalized Sackin's index I_S^*.
# ==============================================================================

"""
    solve_forward_epidemic(zp::ZoonosisParams; tmax::Float64=35.0, dt::Float64=0.01, solver=Tsit5())

Solve the forward epidemiological ODEs on `[0, tmax]` using `OrdinaryDiffEqTsit5`
and return `(sol = sol, df = df, t_peak = t_peak)`.
"""
function solve_forward_epidemic(
    zp::ZoonosisParams;
    tmax::Float64 = 35.0,
    dt::Float64 = 0.01,
    solver = Tsit5()
)
    function fwd_rhs!(du, u, p, t)
        X1, Y1, X2, Y2, C11, C12, C21, C22 = u
        f11 = p.beta * X1 * Y1
        f12 = p.epsilon * p.beta * X2 * Y1
        f21 = p.psi * p.beta * X1 * Y2
        f22 = p.phi * p.beta * X2 * Y2
        du[1] = p.nu1 - f11 - f21 - p.mu1 * X1
        du[2] = f11 + f21 - (p.gamma1 + p.mu1 + p.m1) * Y1
        du[3] = p.nu2 - f12 - f22 - p.mu2 * X2
        du[4] = f12 + f22 - (p.gamma2 + p.mu2 + p.m2) * Y2
        du[5] = f11
        du[6] = f12
        du[7] = f21
        du[8] = f22
        return nothing
    end

    u0 = [zp.X1_0, zp.Y1_0, zp.X2_0, zp.Y2_0, 0.0, 0.0, 0.0, 0.0]
    saveat = collect(0.0:dt:tmax)
    prob = ODEProblem(fwd_rhs!, u0, (0.0, tmax), zp)
    sol = solve(prob, solver; saveat=saveat, reltol=1e-10, abstol=1e-12)

    y1_vec = [u[2] for u in sol.u]
    idx_peak = argmax(y1_vec)
    t_peak = sol.t[idx_peak]

    df = DataFrame(
        time = sol.t,
        X1 = [u[1] for u in sol.u],
        Y1 = [u[2] for u in sol.u],
        X2 = [u[3] for u in sol.u],
        Y2 = [u[4] for u in sol.u],
        cum_f11 = [u[5] for u in sol.u],
        cum_f12 = [u[6] for u in sol.u],
        cum_f21 = [u[7] for u in sol.u],
        cum_f22 = [u[8] for u in sol.u],
        f11 = [zp.beta * u[1] * u[2] for u in sol.u],
        f12 = [zp.epsilon * zp.beta * u[3] * u[2] for u in sol.u],
        f21 = [zp.psi * zp.beta * u[1] * u[4] for u in sol.u],
        f22 = [zp.phi * zp.beta * u[3] * u[4] for u in sol.u]
    )
    return (sol = sol, df = df, t_peak = t_peak)
end

@inline function _general_parsimony_4state_derivs(
    A11::Float64, A12::Float64, A21::Float64, A22::Float64,
    I1::Float64, I2::Float64,
    f11::Float64, f12::Float64, f21::Float64, f22::Float64,
    w::Float64
)
    r11 = A11 / I1
    r12 = A12 / I1
    r21 = A21 / I2
    r22 = A22 / I2

    p1, p2 = r11 + r12, r21 + r22
    migration1 = f12 * r21 * (1 - p1) - f21 * r11 * (1 - p2)
    migration2 = f12 * r22 * (1 - p1) - f21 * r12 * (1 - p2)
    switch1 = f12 * ((1 - w) * r12 * r21 - w * r11 * r22)
    switch2 = f21 * ((1 - w) * r22 * r11 - w * r21 * r12)
    dA11 = -f11 * r11^2 - 2w * f11 * r11 * r12 +
        switch1 - f21 * r11 * p2 + migration1
    dA12 = -f11 * r12^2 - 2(1 - w) * f11 * r11 * r12 -
        switch1 - f21 * r12 * p2 + migration2
    dA21 = -f22 * r21^2 - 2w * f22 * r21 * r22 +
        switch2 - f12 * r21 * p1 - migration1
    dA22 = -f22 * r22^2 - 2(1 - w) * f22 * r21 * r22 -
        switch2 - f12 * r22 * p1 - migration2
    dJ = 2f11 * r11 * r12 + 2f22 * r21 * r22 +
        (f12 + f21) * (r11 * r22 + r12 * r21)
    return (dA11, dA12, dA21, dA22, dJ)
end

@inline function _fitch_parsimony_6state_derivs(
    F11::Float64, F12::Float64, F1s::Float64,
    F21::Float64, F22::Float64, F2s::Float64,
    I1::Float64, I2::Float64,
    f11::Float64, f12::Float64, f21::Float64, f22::Float64
)
    p11 = F11 / I1; p12 = F12 / I1; p1s = F1s / I1; p1 = p11 + p12 + p1s
    p21 = F21 / I2; p22 = F22 / I2; p2s = F2s / I2; p2 = p21 + p22 + p2s

    dF11 = -f11 * p11 * (p11 + 2.0 * p12) +
           f12 * p21 * (1.0 - p1) - f21 * p11 * (1.0 - p2) +
           f12 * (p1s * p21 - p11 * p22) -
           f21 * p11 * p2

    dF12 = -f11 * p12 * (p12 + 2.0 * p11) +
           f12 * p22 * (1.0 - p1) - f21 * p12 * (1.0 - p2) +
           f12 * (p1s * p22 - p12 * p21) -
           f21 * p12 * p2

    dF1s = f11 * (2.0 * p11 * p12 - 2.0 * p1s * p1 + p1s^2) +
           f12 * p2s * (1.0 - p1) - f21 * p1s * (1.0 - p2) +
           f12 * (p11 * p22 + p12 * p21 - p1s * p21 - p1s * p22) -
           f21 * p1s * p2

    dF21 = -f22 * p21 * (p21 + 2.0 * p22) -
           f12 * p21 * (1.0 - p1) + f21 * p11 * (1.0 - p2) -
           f12 * p21 * p1 +
           f21 * (p2s * p11 - p21 * p12)

    dF22 = -f22 * p22 * (p22 + 2.0 * p21) -
           f12 * p22 * (1.0 - p1) + f21 * p12 * (1.0 - p2) -
           f12 * p22 * p1 +
           f21 * (p2s * p12 - p22 * p11)

    dF2s = f22 * (2.0 * p21 * p22 - 2.0 * p2s * p2 + p2s^2) -
           f12 * p2s * (1.0 - p1) + f21 * p1s * (1.0 - p2) -
           f12 * p2s * p1 +
           f21 * (p21 * p12 + p22 * p11 - p2s * p11 - p2s * p12)

    dJF = 2.0 * f11 * p11 * p12 + 2.0 * f22 * p21 * p22 +
          (f12 + f21) * (p11 * p22 + p12 * p21)

    return (dF11, dF12, dF1s, dF21, dF22, dF2s, dJF)
end

function _zoonosis_coalescent_rhs!(du, y, p, uf;
    positive = x -> max(x, 0.0),
    descendant_mean = (x, a) -> a > 1e-9 ? x / a : 0.0)
        X1_f = positive(uf[1])
        I1   = uf[2]
        X2_f = positive(uf[3])
        I2   = uf[4]

        f11 = positive(p.beta * X1_f * I1)
        f12 = positive(p.epsilon * p.beta * X2_f * I1)
        f21 = positive(p.psi * p.beta * X1_f * I2)
        f22 = positive(p.phi * p.beta * X2_f * I2)

        a1 = positive(y[1])
        a2 = positive(y[2])
        p1 = a1 / I1
        p2 = a2 / I2

        L11 = positive(y[3])
        L21 = positive(y[4])
        L12 = positive(y[5])
        L22 = positive(y[6])

        X1_s = positive(y[10])
        X2_s = positive(y[11])
        x1_ratio = descendant_mean(X1_s, a1)
        x2_ratio = descendant_mean(X2_s, a2)

        # 1. True deme lineage counts (Volz 2012 Eqs. 38-39)
        dA1 = -f11 * p1^2 - f21 * p1 * p2 + f12 * p2 * (1.0 - p1) - f21 * p1 * (1.0 - p2)
        dA2 = -f22 * p2^2 - f12 * p1 * p2 + f21 * p1 * (1.0 - p2) - f12 * p2 * (1.0 - p1)

        # 2. External leaf lineages L_{kl} and cherries C11, C12, C22
        dL11 = -(L11 / I1) * (2.0 * f11 * p1 + (f12 + f21) * p2) - (L11 / I1) * f21 * (1.0 - p2) + (L21 / I2) * f12 * (1.0 - p1)
        dL21 = -(L21 / I2) * (2.0 * f22 * p2 + (f12 + f21) * p1) - (L21 / I2) * f12 * (1.0 - p1) + (L11 / I1) * f21 * (1.0 - p2)
        dL12 = -(L12 / I1) * (2.0 * f11 * p1 + (f12 + f21) * p2) - (L12 / I1) * f21 * (1.0 - p2) + (L22 / I2) * f12 * (1.0 - p1)
        dL22 = -(L22 / I2) * (2.0 * f22 * p2 + (f12 + f21) * p1) - (L22 / I2) * f12 * (1.0 - p1) + (L12 / I1) * f21 * (1.0 - p2)

        dC11 = f11 * (L11 / I1)^2 + (f12 + f21) * (L11 / I1) * (L21 / I2) + f22 * (L21 / I2)^2
        dC12 = 2.0 * f11 * (L11 / I1) * (L12 / I1) + (f12 + f21) * ((L11 / I1) * (L22 / I2) + (L12 / I1) * (L21 / I2)) + 2.0 * f22 * (L21 / I2) * (L22 / I2)
        dC22 = f22 * (L22 / I2)^2 + (f12 + f21) * (L12 / I1) * (L22 / I2) + f11 * (L12 / I1)^2

        # 3. Sackin's index auxiliary variables (X1_s, X2_s, K, K1, K2)
        dX1_s = f12 * p2 * x2_ratio - f21 * p1 * x1_ratio
        dX2_s = f21 * p1 * x1_ratio - f12 * p2 * x2_ratio
        dK    = 2.0 * f11 * p1^2 * x1_ratio + 2.0 * f22 * p2^2 * x2_ratio + (f12 + f21) * p1 * p2 * (x1_ratio + x2_ratio)
        dK1   = 2.0 * f11 * p1^2 * x1_ratio
        dK2   = (f12 + f21) * p1 * p2 * (x1_ratio + x2_ratio) + 2.0 * f22 * p2^2 * x2_ratio

        # 4. Prior 1: Unidirectional Animal -> Human Camin-Sokal (A11, A12, A21, A22; J_p1)
        p11 = positive(y[15]) / I1; p12 = positive(y[16]) / I1
        p21 = positive(y[17]) / I2; p22 = positive(y[18]) / I2
        dA11 = -f11 * p11^2 + f12 * p12 * p21 - f21 * p11 * (p21 + p22) + f12 * p21 * (1.0 - p1) - f21 * p11 * (1.0 - p2)
        dA12 = -f11 * p12^2 - 2.0 * f11 * p11 * p12 - f12 * p12 * p21 - f21 * p12 * (p21 + p22) + f12 * p22 * (1.0 - p1) - f21 * p12 * (1.0 - p2)
        dA21 = -f22 * p21^2 + f21 * p22 * p11 - f12 * p21 * (p11 + p12) + f21 * p11 * (1.0 - p2) - f12 * p21 * (1.0 - p1)
        dA22 = -f22 * p22^2 - 2.0 * f22 * p21 * p22 - f21 * p22 * p11 - f12 * p22 * (p11 + p12) + f21 * p12 * (1.0 - p2) - f12 * p22 * (1.0 - p1)
        dJ_p1 = 2.0 * f11 * p11 * p12 + 2.0 * f22 * p21 * p22 + (f12 + f21) * (p11 * p22 + p12 * p21)

        # 5. Prior 2: Unidirectional Human -> Animal Reverse Camin-Sokal (B11, B12, B21, B22; J_p2)
        q11 = positive(y[19]) / I1; q12 = positive(y[20]) / I1
        q21 = positive(y[21]) / I2; q22 = positive(y[22]) / I2
        dB11 = -f11 * q11^2 - 2.0 * f11 * q11 * q12 - f12 * q11 * q22 - f21 * q11 * (q21 + q22) + f12 * q21 * (1.0 - p1) - f21 * q11 * (1.0 - p2)
        dB12 = -f11 * q12^2 + f12 * q11 * q22 - f21 * q12 * (q21 + q22) + f12 * q22 * (1.0 - p1) - f21 * q12 * (1.0 - p2)
        dB21 = -f22 * q21^2 - 2.0 * f22 * q21 * q22 - f21 * q21 * q12 - f12 * q21 * (q11 + q12) + f21 * q11 * (1.0 - p2) - f12 * q21 * (1.0 - p1)
        dB22 = -f22 * q22^2 + f21 * q21 * q12 - f12 * q22 * (q11 + q12) + f21 * q12 * (1.0 - p2) - f12 * q22 * (1.0 - p1)
        dJ_p2 = 2.0 * f11 * q11 * q12 + 2.0 * f22 * q21 * q22 + (f12 + f21) * (q11 * q22 + q12 * q21)

        # 6. Prior 3: Bidirectional Fitch Parsimony (F11, F12, F1s, F21, F22, F2s; J_fitch)
        r11 = positive(y[23]) / I1; r12 = positive(y[24]) / I1; r1s = positive(y[25]) / I1
        r21 = positive(y[26]) / I2; r22 = positive(y[27]) / I2; r2s = positive(y[28]) / I2
        r1 = r11 + r12 + r1s; r2 = r21 + r22 + r2s
        dF11 = f11 * (r11^2 + 2.0*r11*r1s - r11*r1) + f12 * (r11*r21 + r11*r2s + r1s*r21 - r11*r2) - f11*r1*r11 - f21*r2*r11 + f12*r21*(1.0 - r1) - f21*r11*(1.0 - r2)
        dF12 = f11 * (r12^2 + 2.0*r12*r1s - r12*r1) + f12 * (r12*r22 + r12*r2s + r1s*r22 - r12*r2) - f11*r1*r12 - f21*r2*r12 + f12*r22*(1.0 - r1) - f21*r12*(1.0 - r2)
        dF1s = f11 * (r1s^2 + 2.0*r11*r12 - r1s*r1) + f12 * (r1s*r2s + r11*r22 + r12*r21 - r1s*r2) - f11*r1*r1s - f21*r2*r1s + f12*r2s*(1.0 - r1) - f21*r1s*(1.0 - r2)
        dF21 = f22 * (r21^2 + 2.0*r21*r2s - r21*r2) + f21 * (r21*r11 + r21*r1s + r2s*r11 - r21*r1) - f22*r2*r21 - f12*r1*r21 + f21*r11*(1.0 - r2) - f12*r21*(1.0 - r1)
        dF22 = f22 * (r22^2 + 2.0*r22*r2s - r22*r2) + f21 * (r22*r12 + r22*r1s + r2s*r12 - r22*r1) - f22*r2*r22 - f12*r1*r22 + f21*r12*(1.0 - r2) - f12*r22*(1.0 - r1)
        dF2s = f22 * (r2s^2 + 2.0*r21*r22 - r2s*r2) + f21 * (r2s*r1s + r21*r12 + r22*r11 - r2s*r1) - f22*r2*r2s - f12*r1*r2s + f21*r1s*(1.0 - r2) - f12*r2s*(1.0 - r1)
        dJ_fitch = 2.0 * f11 * r11 * r12 + 2.0 * f22 * r21 * r22 + (f12 + f21) * (r11 * r22 + r12 * r21)

        dJs   = f12 * p2
        dJs21 = f21 * p1
        dTL    = a1 + a2
        dTLA1  = a1
        dTLA11 = positive(y[15]) + positive(y[17])
        dTLA2  = a2

        du[1] = dA1; du[2] = dA2
        du[3] = dL11; du[4] = dL21; du[5] = dL12; du[6] = dL22
        du[7] = dC11; du[8] = dC12; du[9] = dC22
        du[10] = dX1_s; du[11] = dX2_s; du[12] = dK; du[13] = dK1; du[14] = dK2
        du[15] = dA11; du[16] = dA12; du[17] = dA21; du[18] = dA22
        du[19] = dB11; du[20] = dB12; du[21] = dB21; du[22] = dB22
        du[23] = dF11; du[24] = dF12; du[25] = dF1s; du[26] = dF21; du[27] = dF22; du[28] = dF2s
        du[29] = dJs; du[30] = dJs21; du[31] = dJ_p1; du[32] = dJ_p2; du[33] = dJ_fitch
        du[34] = dTL; du[35] = dTLA1; du[36] = dTLA11; du[37] = dTLA2
        return nothing
end

"""
    solve_coalescent_odes(zp::ZoonosisParams; n1::Float64=100.0, n2::Float64=100.0, T_sample::Float64=zp.t_peak, dt::Float64=0.01, solver=Tsit5())

Solve the two-host, 37-state Frost/Volz backward coalescent ODE system from `s = 0`
(at `t = T_sample`) to `s = T_sample` (at `t = 0`), returning
`(summary = ..., traj = DataFrame, fwd = fwd)`. For selectable multi-deme equations,
use the `TFGY` or `ReactionSystem` overload.
"""
function solve_coalescent_odes(
    zp::ZoonosisParams;
    n1::Float64 = 100.0,
    n2::Float64 = 100.0,
    T_sample::Float64 = zp.t_peak,
    dt::Float64 = 0.01,
    solver = Tsit5()
)
    fwd = solve_forward_epidemic(zp; tmax = max(T_sample + 1.0, 35.0), dt = dt, solver = solver)
    fwd_sol = fwd.sol
    u_T = fwd_sol(T_sample)
    Y1_T = max(u_T[2], 1e-9)
    Y2_T = max(u_T[4], 1e-9)
    n1_eff = n1
    n2_eff = n2
    n_total = n1_eff + n2_eff
    function bwd_rhs!(du, y, p, s)
        uf = fwd_sol(max(0.0, T_sample - s))
        uf_eff = [max(uf[1], 1e-12), max(uf[2], 1e-9),
            max(uf[3], 1e-12), max(uf[4], 1e-9)]
        return _zoonosis_coalescent_rhs!(du, y, p, uf_eff)
    end

    y0 = zeros(Float64, 37)
    y0[1] = n1_eff; y0[2] = n2_eff
    y0[3] = n1_eff; y0[6] = n2_eff
    y0[10] = n1_eff; y0[11] = n2_eff
    y0[15] = n1_eff; y0[18] = n2_eff
    y0[19] = n1_eff; y0[22] = n2_eff
    y0[23] = n1_eff; y0[27] = n2_eff

    s_grid = collect(0.0:dt:T_sample)
    if s_grid[end] < T_sample
        push!(s_grid, T_sample)
    end
    prob = ODEProblem(bwd_rhs!, y0, (0.0, T_sample), zp)
    sol = solve(prob, solver; saveat=s_grid, reltol=1e-9, abstol=1e-11)

    Atot = [u[1] + u[2] for u in sol.u]
    idx_tmrca = argmin(abs.(Atot .- 1.0))
    yf = sol.u[idx_tmrca]

    C11 = yf[7]; C12 = yf[8]; C22 = yf[9]
    denom_r = (C11 + 0.5 * C12) * (C22 + 0.5 * C12)
    r_assort = denom_r > 0.0 ? (C11 * C22 - 0.25 * C12^2) / denom_r : 0.0

    IS_raw = yf[12]
    n_int = max(2, round(Int, n_total))
    E_IS_yule = 2.0 * n_total * sum(1.0 / j for j in 2:n_int)
    IS_norm = (IS_raw - E_IS_yule) / E_IS_yule

    J_true_12    = yf[29]
    J_true_21    = yf[30]
    J_true_total = J_true_12 + J_true_21
    J_pars_1to2  = yf[31]
    J_pars_2to1  = yf[32]
    J_pars_fitch = yf[33]

    summary = (
        T_sample = T_sample,
        n1 = n1_eff,
        n2 = n2_eff,
        Y1_T = Y1_T,
        Y2_T = Y2_T,
        J_true_12 = J_true_12,
        J_true_21 = J_true_21,
        J_true_total = J_true_total,
        J_true = J_true_12,
        J_pars = J_pars_1to2,
        J_pars_1to2 = J_pars_1to2,
        J_pars_2to1 = J_pars_2to1,
        J_pars_fitch = J_pars_fitch,
        ratio_pars_true = J_true_12 > 0 ? J_pars_1to2 / J_true_12 : NaN,
        ratio_p1_total = J_true_total > 0 ? J_pars_1to2 / J_true_total : NaN,
        ratio_p2_total = J_true_total > 0 ? J_pars_2to1 / J_true_total : NaN,
        ratio_fitch_total = J_true_total > 0 ? J_pars_fitch / J_true_total : NaN,
        C11 = C11,
        C12 = C12,
        C22 = C22,
        r_assort = r_assort,
        IS_raw = IS_raw,
        IS_yule = E_IS_yule,
        IS_norm = IS_norm,
        L11_pure = yf[36],
        L1_true = yf[35],
        L2_true = yf[37],
        L_total = yf[34],
        prop_human_branch = yf[37] / max(yf[34], 1e-12)
    )

    traj = DataFrame(
        s = sol.t,
        t_fwd = T_sample .- sol.t,
        A1 = [u[1] for u in sol.u],
        A2 = [u[2] for u in sol.u],
        C11 = [u[7] for u in sol.u],
        C12 = [u[8] for u in sol.u],
        C22 = [u[9] for u in sol.u],
        IS_raw = [u[12] for u in sol.u],
        A11 = [u[15] for u in sol.u],
        A12 = [u[16] for u in sol.u],
        A21 = [u[17] for u in sol.u],
        A22 = [u[18] for u in sol.u],
        B11 = [u[19] for u in sol.u],
        B12 = [u[20] for u in sol.u],
        B21 = [u[21] for u in sol.u],
        B22 = [u[22] for u in sol.u],
        F11 = [u[23] for u in sol.u],
        F12 = [u[24] for u in sol.u],
        F1s = [u[25] for u in sol.u],
        F21 = [u[26] for u in sol.u],
        F22 = [u[27] for u in sol.u],
        F2s = [u[28] for u in sol.u],
        J_true_12 = [u[29] for u in sol.u],
        J_true_21 = [u[30] for u in sol.u],
        J_true_total = [u[29] + u[30] for u in sol.u],
        J_true = [u[29] for u in sol.u],
        J_pars = [u[31] for u in sol.u],
        J_pars_1to2 = [u[31] for u in sol.u],
        J_pars_2to1 = [u[32] for u in sol.u],
        J_pars_fitch = [u[33] for u in sol.u],
        L1_true = [u[35] for u in sol.u],
        L2_true = [u[37] for u in sol.u]
    )

    return (summary = summary, traj = traj, fwd = fwd)
end

"""
    compute_canonical_summary(; n1=100.0, n2=100.0, T_sample=12.85) -> DataFrame

Compute the analytical coalescent ODE summary table across Models 1, 2, and 3.
"""
function compute_canonical_summary(; n1::Float64=100.0, n2::Float64=100.0, T_sample::Float64=12.85)
    models = [
        ("Model 1 (Pure Spillover)", model1_params(t_peak=T_sample)),
        ("Model 2 (With Human-to-Human Spread)", model2_params(t_peak=T_sample)),
        ("Model 3 (With Reverse Zoonosis)", model3_params(t_peak=T_sample))
    ]
    rows = NamedTuple[]
    for (mname, zp) in models
        res = solve_coalescent_odes(zp; n1=n1, n2=n2, T_sample=T_sample)
        push!(rows, merge((model = mname, epsilon = zp.epsilon, phi = zp.phi, psi = zp.psi), res.summary))
    end
    return DataFrame(rows)
end

"""
    sweep_sampling_time(; T_vals=4.0:0.5:25.0, n1=100.0, n2=100.0) -> DataFrame

Sweep sampling time `T` across Models 1, 2, and 3.
"""
function sweep_sampling_time(; T_vals = 4.0:0.5:25.0, n1::Float64 = 100.0, n2::Float64 = 100.0)
    models = [
        ("Model 1", model1_params()),
        ("Model 2", model2_params()),
        ("Model 3", model3_params())
    ]
    rows = NamedTuple[]
    for (mname, zp) in models
        for T in T_vals
            res = solve_coalescent_odes(zp; n1=n1, n2=n2, T_sample=Float64(T), dt=0.02)
            push!(rows, merge((model = mname,), res.summary))
        end
    end
    return DataFrame(rows)
end

"""
    sweep_sampling_fraction(; n_vals=[10, 20, 40, 60, 80, 100, 140, 180, 220, 260, 300], T_sample=12.85) -> (df_n2, df_n1)

Sweep human sample size `n2` (fixing `n1=100`) and reservoir sample size `n1` (fixing `n2=100`)
across Models 1, 2, and 3.
"""
function sweep_sampling_fraction(;
    n_vals = [10.0, 20.0, 40.0, 60.0, 80.0, 100.0, 140.0, 180.0, 220.0, 260.0, 300.0],
    T_sample::Float64 = 12.85
)
    models = [
        ("Model 1", model1_params()),
        ("Model 2", model2_params()),
        ("Model 3", model3_params())
    ]
    rows_n2 = NamedTuple[]
    rows_n1 = NamedTuple[]
    for (mname, zp) in models
        for nv in n_vals
            r2 = solve_coalescent_odes(zp; n1=100.0, n2=Float64(nv), T_sample=T_sample, dt=0.02)
            push!(rows_n2, merge((model = mname,), r2.summary))
            r1 = solve_coalescent_odes(zp; n1=Float64(nv), n2=100.0, T_sample=T_sample, dt=0.02)
            push!(rows_n1, merge((model = mname,), r1.summary))
        end
    end
    return (df_n2 = DataFrame(rows_n2), df_n1 = DataFrame(rows_n1))
end

"""
    sweep_R22(; R22_vals=0.0:0.25:5.0, T_sample=12.85) -> DataFrame

Sweep human-to-human reproduction number `R22` for `epsilon in (0.01, 0.1)`.
"""
function sweep_R22(; R22_vals = 0.0:0.25:5.0, T_sample::Float64 = 12.85)
    base = model1_params()
    denom2 = base.gamma2 + base.mu2 + base.m2
    rows = NamedTuple[]
    for eps_val in (0.01, 0.1)
        for r22 in R22_vals
            phi_val = r22 * denom2 / (base.beta * base.N2_0)
            zp = ZoonosisParams(epsilon=eps_val, phi=phi_val, psi=0.0, t_peak=T_sample)
            res = solve_coalescent_odes(zp; n1=100.0, n2=100.0, T_sample=T_sample, dt=0.02)
            push!(rows, merge((R22 = Float64(r22), epsilon = eps_val, phi = phi_val), res.summary))
        end
    end
    return DataFrame(rows)
end

"""
    sweep_psi(; psi_vals=[0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0], T_sample=12.85) -> DataFrame

Sweep reverse zoonosis multiplier `psi` across `[0, 8]`.
"""
function sweep_psi(;
    psi_vals = [0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0],
    T_sample::Float64 = 12.85
)
    rows = NamedTuple[]
    for psi_val in psi_vals
        zp = model2_params(psi=Float64(psi_val), t_peak=T_sample)
        res = solve_coalescent_odes(zp; n1=100.0, n2=100.0, T_sample=T_sample, dt=0.02)
        push!(rows, merge((psi = Float64(psi_val), epsilon = zp.epsilon, phi = zp.phi), res.summary))
    end
    return DataFrame(rows)
end
