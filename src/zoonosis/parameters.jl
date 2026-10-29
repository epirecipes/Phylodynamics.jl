"""
    ZoonosisParams

Epidemiological and sampling parameters for the two-host SIRS zoonosis model
(reservoir host deme 1: `X1, Y1`; human host deme 2: `X2, Y2`).
"""
Base.@kwdef struct ZoonosisParams
    N1_0::Float64   = 10000.0
    X1_0::Float64   = 9999.0
    Y1_0::Float64   = 1.0
    N2_0::Float64   = 1000.0
    X2_0::Float64   = 1000.0
    Y2_0::Float64   = 0.0
    nu1::Float64    = 2.74
    mu1::Float64    = 2.74e-4
    gamma1::Float64 = 0.1
    m1::Float64     = 0.0
    beta::Float64   = 10.0 * (2.74e-4 + 0.1) / 10000.0
    nu2::Float64    = 1000.0 / (70.0 * 365.0)
    mu2::Float64    = 1.0 / (70.0 * 365.0)
    gamma2::Float64 = 0.01
    m2::Float64     = 0.1
    epsilon::Float64 = 0.1
    phi::Float64     = 0.0
    psi::Float64     = 0.0
    t_peak::Float64  = 12.85
end

"""Canonical Model 1 (Pure Zoonotic Spillover: epsilon=0.1, phi=0, psi=0)."""
model1_params(; kwargs...) = ZoonosisParams(; epsilon=0.1, phi=0.0, psi=0.0, kwargs...)

"""Canonical Model 2 (Sustained Human-to-Human Spread: epsilon=0.01, phi=4.0, psi=0)."""
model2_params(; kwargs...) = ZoonosisParams(; epsilon=0.01, phi=4.0, psi=0.0, kwargs...)

"""Canonical Model 3 (Bidirectional / Reverse Zoonosis: epsilon=0.01, phi=4.0, psi=4.0)."""
model3_params(; kwargs...) = ZoonosisParams(; epsilon=0.01, phi=4.0, psi=4.0, kwargs...)
