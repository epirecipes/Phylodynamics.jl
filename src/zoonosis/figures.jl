# ==============================================================================
# figures.jl
# Generates all 9 summary figures (fig1.pdf - fig9.pdf + PNG) using CairoMakie.jl.
# ==============================================================================

const COL_RES    = Makie.wong_colors()[1] # Blue (Reservoir / Model 1 / Prior 1)
const COL_HUM    = Makie.wong_colors()[6] # Red/Vermilion (Human / Model 2 / Prior 2)
const COL_REV    = Makie.wong_colors()[3] # Green (Model 3 / Prior 3 Fitch)
const COL_BLACK  = RGBf(0.15, 0.15, 0.15)
const COL_PURPLE = Makie.wong_colors()[4]

function _save_both(fig::Figure, outdir::AbstractString, basename::AbstractString)
    mkpath(outdir)
    pdf_path = joinpath(outdir, "$(basename).pdf")
    png_path = joinpath(outdir, "$(basename).png")
    save(pdf_path, fig)
    save(png_path, fig; px_per_unit=2)
    return pdf_path
end

"""
    plot_fig1_epidemic_and_coalescent(res1, res2, res3; outdir="julia_output/figures")

Figure 1: Forward epidemic trajectories (`Y1(t), Y2(t)`) and backward coalescent
parsimony vs true jump trajectories across Models 1, 2, and 3.
"""
function plot_fig1_epidemic_and_coalescent(res1, res2, res3; outdir::AbstractString="julia_output/figures")
    fig = Figure(size=(1200, 720), fontsize=13)
    models = [("Model 1: Pure Spillover (φ=0, ψ=0)", res1),
              ("Model 2: Human Spread (φ=4, ψ=0)", res2),
              ("Model 3: Reverse Zoonosis (φ=4, ψ=4)", res3)]

    for (col, (title_str, res)) in enumerate(models)
        df_fwd = filter(r -> r.time <= 30.0, res.fwd.df)
        ax1 = Axis(fig[1, col];
            title = "($(Char('A' + col - 1))) $title_str",
            xlabel = "Forward time t",
            ylabel = col == 1 ? "Infected prevalence Y_k(t)" : "",
            yscale = log10
        )
        lines!(ax1, df_fwd.time, max.(df_fwd.Y1, 1e-2); color=COL_RES, linewidth=2.5, label="Reservoir Y₁(t)")
        lines!(ax1, df_fwd.time, max.(df_fwd.Y2, 1e-2); color=COL_HUM, linewidth=2.5, label="Human Y₂(t)")
        vlines!(ax1, [res.summary.T_sample]; color=:gray40, linestyle=:dash, linewidth=1.5)
        ylims!(ax1, 1e-2, 2e4)
        if col == 1
            axislegend(ax1; position=:rb, framevisible=false, labelsize=11)
        end

        tr = res.traj
        ax2 = Axis(fig[2, col];
            title = "($(Char('D' + col - 1))) Coalescent Jumps (Backward from T=12.85)",
            xlabel = "Forward time t = T - s",
            ylabel = col == 1 ? "Cumulative jumps in genealogy" : ""
        )
        lines!(ax2, tr.t_fwd, tr.J_true_12; color=COL_BLACK, linewidth=2.5, label="True 1→2 (J*₁₂)")
        if maximum(tr.J_true_21) > 0.01
            lines!(ax2, tr.t_fwd, tr.J_true_21; color=COL_PURPLE, linewidth=2.2, linestyle=:dot, label="True 2→1 (J*₂₁)")
        end
        lines!(ax2, tr.t_fwd, tr.J_pars_1to2; color=COL_RES, linewidth=2.2, linestyle=:dash, label="Prior 1 (1→2 CS)")
        lines!(ax2, tr.t_fwd, tr.J_pars_2to1; color=COL_HUM, linewidth=2.2, linestyle=:dashdot, label="Prior 2 (2→1 CS)")
        lines!(ax2, tr.t_fwd, tr.J_pars_fitch; color=COL_REV, linewidth=2.2, label="Prior 3 (Fitch)")
        if col == 1
            axislegend(ax2; position=:lt, framevisible=false, labelsize=10)
        end
    end

    _save_both(fig, outdir, "fig1")
    return fig
end

"""
    plot_fig2_simulated_trees(t1::PhyloTree, t2::PhyloTree, t3::PhyloTree; outdir="julia_output/figures")

Figure 2: Simulated phylogenies across Models 1, 2, and 3 colored by host state
(Blue = Reservoir `Y1`, Red = Human `Y2`).
"""
function plot_fig2_simulated_trees(t1::PhyloTree, t2::PhyloTree, t3::PhyloTree; outdir::AbstractString="julia_output/figures")
    fig = Figure(size=(1250, 520), fontsize=13)
    trees = [("Model 1: Pure Spillover\n(Interspersed Singleton Human Tips)", t1),
             ("Model 2: With Human Spread\n(Monophyletic Human Clades)", t2),
             ("Model 3: With Reverse Zoonosis\n(Bidirectional Host Re-entry)", t3)]

    for (col, (title_str, tr)) in enumerate(trees)
        bt, x, y = _tree_layout_coords(tr)
        pars = parsimony_counts(bt)
        sub_str = @sprintf("J*₁₂=%d, J*₂₁=%d | J⁽⁰⁾=%d, J⁽¹⁾=%d, J_F=%d",
            bt.true_jumps_12, bt.true_jumps_21, pars.J_pars_1to2, pars.J_pars_2to1, pars.J_pars_fitch)
        ax = Axis(fig[1, col];
            title = "($(Char('A' + col - 1))) $title_str\n$sub_str",
            xlabel = "Time from root (years)",
            ylabel = col == 1 ? "Sampled lineages" : "",
            yticklabelsvisible = false,
            yticksvisible = false
        )
        for e in 1:size(bt.edge, 1)
            p = bt.edge[e, 1]
            c = bt.edge[e, 2]
            st_c = bt.node_state[c] > 0 ? bt.node_state[c] : (c <= ntips(bt) ? bt.tip_state[c] : 1)
            col_e = st_c == 2 ? COL_HUM : COL_RES
            lines!(ax, [x[p], x[p], x[c]], [y[p], y[c], y[c]]; color=col_e, linewidth=1.2)
        end
        n = ntips(bt)
        idx1 = findall(==(1), bt.tip_state)
        idx2 = findall(==(2), bt.tip_state)
        !isempty(idx1) && scatter!(ax, x[idx1], y[idx1]; color=COL_RES, markersize=5, label="Reservoir tip (Y₁)")
        !isempty(idx2) && scatter!(ax, x[idx2], y[idx2]; color=COL_HUM, markersize=5, label="Human tip (Y₂)")
        if col == 1
            axislegend(ax; position=:lt, framevisible=true, labelsize=11)
        end
    end

    _save_both(fig, outdir, "fig2")
    return fig
end

"""
    plot_fig3_sampling_time_bias(df_time::DataFrame; outdir="julia_output/figures")

Figure 3: Parsimony jump counts and recovery ratios across sampling times `T in [4, 25]`
for Models 1, 2, and 3.
"""
function plot_fig3_sampling_time_bias(df_time::DataFrame; outdir::AbstractString="julia_output/figures")
    fig = Figure(size=(1200, 440), fontsize=13)
    for (col, mname) in enumerate(["Model 1", "Model 2", "Model 3"])
        sub = filter(r -> r.model == mname, df_time)
        ax = Axis(fig[1, col];
            title = "($(Char('A' + col - 1))) $mname: Jump Counts vs Sampling Time T",
            xlabel = "Sampling time T",
            ylabel = col == 1 ? "Expected jump count" : ""
        )
        lines!(ax, sub.T_sample, sub.J_true_total; color=COL_BLACK, linewidth=2.6, label="True total J*")
        lines!(ax, sub.T_sample, sub.J_pars_1to2; color=COL_RES, linewidth=2.2, linestyle=:dash, label="Prior 1 (1→2 CS)")
        lines!(ax, sub.T_sample, sub.J_pars_2to1; color=COL_HUM, linewidth=2.2, linestyle=:dashdot, label="Prior 2 (2→1 CS)")
        lines!(ax, sub.T_sample, sub.J_pars_fitch; color=COL_REV, linewidth=2.2, label="Prior 3 (Fitch)")
        vlines!(ax, [12.85]; color=:gray50, linestyle=:dot, linewidth=1.5)
        if col == 1
            axislegend(ax; position=:lt, framevisible=false, labelsize=10)
        end
    end
    _save_both(fig, outdir, "fig3")
    return fig
end

"""
    plot_fig4_sampling_fraction_bias(df_n2::DataFrame, df_n1::DataFrame; outdir="julia_output/figures")

Figure 4: Effect of human sampling effort `n2` and reservoir sampling effort `n1`
on parsimony recovery ratios and cherry assortativity `r`.
"""
function plot_fig4_sampling_fraction_bias(df_n2::DataFrame, df_n1::DataFrame; outdir::AbstractString="julia_output/figures")
    fig = Figure(size=(1100, 720), fontsize=13)
    colors = Dict("Model 1" => COL_RES, "Model 2" => COL_HUM, "Model 3" => COL_REV)

    ax1 = Axis(fig[1, 1]; title="(A) Parsimony Ratio J_Fitch / J* vs Human Sample n₂", xlabel="Human sample size n₂ (n₁=100)", ylabel="Recovery ratio J_Fitch / J*")
    ax2 = Axis(fig[1, 2]; title="(B) Parsimony Ratio J_Fitch / J* vs Reservoir Sample n₁", xlabel="Reservoir sample size n₁ (n₂=100)", ylabel="Recovery ratio J_Fitch / J*")
    ax3 = Axis(fig[2, 1]; title="(C) Cherry Assortativity r vs Human Sample n₂", xlabel="Human sample size n₂ (n₁=100)", ylabel="Cherry assortativity r")
    ax4 = Axis(fig[2, 2]; title="(D) Prior 1 vs Prior 2 Asymmetry vs Human Sample n₂", xlabel="Human sample size n₂ (n₁=100)", ylabel="J⁽¹⁾ / J⁽⁰⁾ ratio")

    for mname in ["Model 1", "Model 2", "Model 3"]
        s2 = filter(r -> r.model == mname, df_n2)
        s1 = filter(r -> r.model == mname, df_n1)
        c = colors[mname]
        lines!(ax1, s2.n2, s2.ratio_fitch_total; color=c, linewidth=2.5, label=mname)
        lines!(ax2, s1.n1, s1.ratio_fitch_total; color=c, linewidth=2.5, label=mname)
        lines!(ax3, s2.n2, s2.r_assort; color=c, linewidth=2.5, label=mname)
        lines!(ax4, s2.n2, s2.J_pars_2to1 ./ max.(s2.J_pars_1to2, 1e-6); color=c, linewidth=2.5, label=mname)
    end
    hlines!(ax4, [1.0]; color=:gray50, linestyle=:dash)
    axislegend(ax1; position=:rb, framevisible=false)
    _save_both(fig, outdir, "fig4")
    return fig
end

"""
    plot_fig5_R22_sweep(df_r22::DataFrame; outdir="julia_output/figures")

Figure 5: Effect of human-to-human reproduction number `R22 in [0, 5]` on true
spillover jumps vs parsimony under Priors 1, 2, and 3.
"""
function plot_fig5_R22_sweep(df_r22::DataFrame; outdir::AbstractString="julia_output/figures")
    fig = Figure(size=(1050, 450), fontsize=13)
    for (col, eps_val) in enumerate((0.01, 0.1))
        sub = filter(r -> r.epsilon == eps_val, df_r22)
        ax = Axis(fig[1, col];
            title = "($(Char('A' + col - 1))) Spillover vs Parsimony across R₂₂ (ε = $eps_val)",
            xlabel = "Human-to-human reproduction number R₂₂",
            ylabel = col == 1 ? "Expected lineage jumps" : ""
        )
        lines!(ax, sub.R22, sub.J_true_12; color=COL_BLACK, linewidth=2.6, label="True J*₁₂")
        lines!(ax, sub.R22, sub.J_pars_1to2; color=COL_RES, linewidth=2.2, linestyle=:dash, label="Prior 1 (1→2 CS)")
        lines!(ax, sub.R22, sub.J_pars_2to1; color=COL_HUM, linewidth=2.2, linestyle=:dashdot, label="Prior 2 (2→1 CS)")
        lines!(ax, sub.R22, sub.J_pars_fitch; color=COL_REV, linewidth=2.2, label="Prior 3 (Fitch)")
        vlines!(ax, [1.0]; color=:gray50, linestyle=:dot, linewidth=1.5)
        if col == 1
            axislegend(ax; position=:rt, framevisible=false)
        end
    end
    _save_both(fig, outdir, "fig5")
    return fig
end

"""
    plot_fig6_parsimony_vs_mk2(df_sim_summary::DataFrame, ode_summary::DataFrame; outdir="julia_output/figures")

Figure 6: Comparison of True jumps, Parsimony (Priors 1, 2, 3), Unidirectional Mk2 ML,
and Bidirectional 2-rate Mk2 ML across Models 1, 2, and 3.
"""
function plot_fig6_parsimony_vs_mk2(df_sim_summary::DataFrame, ode_summary::DataFrame; outdir::AbstractString="julia_output/figures")
    fig = Figure(size=(1150, 500), fontsize=13)

    ax1 = Axis(fig[1, 1];
        title = "(A) Analytical Coalescent ODE vs ReMASTER Simulation Means",
        xlabel = "Epidemiological Regime",
        ylabel = "Jump Count (n₁=n₂=100)",
        xticks = (1:3, ["Model 1\n(Pure Spillover)", "Model 2\n(Human Spread)", "Model 3\n(Reverse Zoonosis)"])
    )
    x_pos = 1:3
    w = 0.14
    barplot!(ax1, x_pos .- 2w, ode_summary.J_true_total; width=w, color=COL_BLACK, label="True Total J* (ODE)")
    barplot!(ax1, x_pos .- 1w, ode_summary.J_pars_1to2; width=w, color=COL_RES, label="Prior 1: 1→2 CS (ODE)")
    barplot!(ax1, x_pos .+ 0w, ode_summary.J_pars_2to1; width=w, color=COL_HUM, label="Prior 2: 2→1 CS (ODE)")
    barplot!(ax1, x_pos .+ 1w, ode_summary.J_pars_fitch; width=w, color=COL_REV, label="Prior 3: Fitch (ODE)")
    barplot!(ax1, x_pos .+ 2w, df_sim_summary.J_ml_bidir_tot; width=w, color=COL_PURPLE, label="2-Rate Mk2 ML (Sim)")
    axislegend(ax1; position=:rt, framevisible=false, labelsize=10)

    ax2 = Axis(fig[1, 2];
        title = "(B) Directional Jump Decomposition (1→2 vs 2→1)",
        xlabel = "Epidemiological Regime",
        ylabel = "Mean Directional Jumps",
        xticks = (1:3, ["Model 1", "Model 2", "Model 3"])
    )
    barplot!(ax2, x_pos .- 1.5w, df_sim_summary.J_true_12; width=w, color=COL_RES, label="True 1→2 (J*₁₂)")
    barplot!(ax2, x_pos .- 0.5w, df_sim_summary.J_true_21; width=w, color=COL_HUM, label="True 2→1 (J*₂₁)")
    barplot!(ax2, x_pos .+ 0.5w, df_sim_summary.J_ml_bidir_12; width=w, color=(COL_RES, 0.55), label="Mk2 ML 1→2")
    barplot!(ax2, x_pos .+ 1.5w, df_sim_summary.J_ml_bidir_21; width=w, color=(COL_HUM, 0.55), label="Mk2 ML 2→1")
    axislegend(ax2; position=:rt, framevisible=false, labelsize=10)

    _save_both(fig, outdir, "fig6")
    return fig
end

"""
    plot_fig7_structured_coalescent_profiles(prof; outdir="julia_output/figures")

Figure 7: Structured coalescent profile log-likelihood curves for `epsilon` (Model 1),
`phi` (Model 2), `psi` (Model 3), and 2D `(epsilon, phi)` surface.
"""
function plot_fig7_structured_coalescent_profiles(prof; outdir::AbstractString="julia_output/figures")
    fig = Figure(size=(1100, 740), fontsize=13)

    ax1 = Axis(fig[1, 1]; title="(A) Model 1 (True φ=0): Profile Log-Likelihood over φ", xlabel="Human transmission multiplier φ (true = 0.0)", ylabel="Δ Log-Likelihood")
    lines!(ax1, prof.profile_eps.phi, prof.profile_eps.delta_ll; color=COL_RES, linewidth=2.6)
    scatter!(ax1, prof.profile_eps.phi, prof.profile_eps.delta_ll; color=COL_RES, markersize=8)
    vlines!(ax1, [0.0]; color=COL_BLACK, linestyle=:dash, linewidth=1.8, label="True φ = 0.0")
    hlines!(ax1, [-1.92]; color=:gray50, linestyle=:dot, label="95% CI threshold (-1.92)")
    ylims!(ax1, -8.0, 0.5)
    axislegend(ax1; position=:lb, framevisible=false)

    ax2 = Axis(fig[1, 2]; title="(B) Model 2 (True φ=4): Profile Log-Likelihood over φ", xlabel="Human transmission multiplier φ (true = 4.0)", ylabel="Δ Log-Likelihood")
    lines!(ax2, prof.profile_phi.phi, prof.profile_phi.delta_ll; color=COL_HUM, linewidth=2.6)
    scatter!(ax2, prof.profile_phi.phi, prof.profile_phi.delta_ll; color=COL_HUM, markersize=8)
    vlines!(ax2, [4.0]; color=COL_BLACK, linestyle=:dash, linewidth=1.8, label="True φ = 4.0")
    hlines!(ax2, [-1.92]; color=:gray50, linestyle=:dot)
    ylims!(ax2, -12.0, 0.5)
    axislegend(ax2; position=:lb, framevisible=false)

    ax3 = Axis(fig[2, 1]; title="(C) Profile Log-Likelihood for Reverse Zoonosis ψ", xlabel="Reverse zoonosis multiplier ψ (Model 3 true ψ=4, Model 2 true ψ=0)", ylabel="Δ Log-Likelihood")
    lines!(ax3, prof.profile_psi.psi, prof.profile_psi.delta_ll; color=COL_REV, linewidth=2.6, label="Model 3 (true ψ = 4.0)")
    scatter!(ax3, prof.profile_psi.psi, prof.profile_psi.delta_ll; color=COL_REV, markersize=8)
    lines!(ax3, prof.profile_psi.psi, prof.profile_psi.delta_ll_m2; color=COL_HUM, linewidth=2.2, linestyle=:dash, label="Model 2 (true ψ = 0.0)")
    vlines!(ax3, [4.0]; color=COL_BLACK, linestyle=:dash, linewidth=1.8)
    hlines!(ax3, [-1.92]; color=:gray50, linestyle=:dot)
    ylims!(ax3, -10.0, 0.5)
    axislegend(ax3; position=:lb, framevisible=false)

    ax4 = Axis(fig[2, 2]; title="(D) Model 2 Joint (ε, φ) Log-Likelihood Slices", xlabel="Human transmission multiplier φ", ylabel="Δ Log-Likelihood")
    for (idx, ev) in enumerate(unique(prof.surface_2d.epsilon))
        sub = filter(r -> r.epsilon == ev, prof.surface_2d)
        lines!(ax4, sub.phi, sub.delta_ll; linewidth=2.0, label="ε = $ev")
    end
    vlines!(ax4, [4.0]; color=COL_BLACK, linestyle=:dash)
    ylims!(ax4, -15.0, 0.5)
    axislegend(ax4; position=:lb, framevisible=false, labelsize=10)

    _save_both(fig, outdir, "fig7")
    return fig
end

"""
    plot_fig8_topology_diagnostics(ode_summary::DataFrame, sim_summary::DataFrame; outdir="julia_output/figures")

Figure 8: Tree topology diagnostics (cherry assortativity `r` and normalized Sackin's index `I_S^*`)
comparing analytical coalescent ODE predictions vs simulated trees across Models 1, 2, and 3.
"""
function plot_fig8_topology_diagnostics(ode_summary::DataFrame, sim_summary::DataFrame; outdir::AbstractString="julia_output/figures")
    fig = Figure(size=(1000, 440), fontsize=13)
    x_pos = 1:3
    w = 0.28

    ax1 = Axis(fig[1, 1];
        title = "(A) Cherry Assortativity Coefficient r",
        xlabel = "Epidemiological Regime",
        ylabel = "Cherry assortativity r",
        xticks = (1:3, ["Model 1", "Model 2", "Model 3"])
    )
    barplot!(ax1, x_pos .- 0.5w, ode_summary.r_assort; width=w, color=COL_RES, label="Coalescent ODE")
    barplot!(ax1, x_pos .+ 0.5w, sim_summary.r_assort; width=w, color=COL_HUM, label="Simulated Trees")
    hlines!(ax1, [0.0]; color=:gray50, linestyle=:dash)
    axislegend(ax1; position=:lt, framevisible=false)

    ax2 = Axis(fig[1, 2];
        title = "(B) Normalized Sackin's Tree Imbalance I_S*",
        xlabel = "Epidemiological Regime",
        ylabel = "Normalized Sackin's index I_S*",
        xticks = (1:3, ["Model 1", "Model 2", "Model 3"])
    )
    barplot!(ax2, x_pos .- 0.5w, ode_summary.IS_norm; width=w, color=COL_RES, label="Coalescent ODE")
    barplot!(ax2, x_pos .+ 0.5w, sim_summary.IS_norm; width=w, color=COL_HUM, label="Simulated Trees")
    hlines!(ax2, [0.0]; color=:gray50, linestyle=:dash)
    axislegend(ax2; position=:lt, framevisible=false)

    _save_both(fig, outdir, "fig8")
    return fig
end

"""
    plot_fig9_reverse_zoonosis_sweep(df_psi::DataFrame; outdir="julia_output/figures")

Figure 9: Impact of reverse zoonosis rate `psi in [0, 8]` on forward (`J*12`) and
reverse (`J*21`) true jumps and parsimony under Priors 1, 2, and 3.
"""
function plot_fig9_reverse_zoonosis_sweep(df_psi::DataFrame; outdir::AbstractString="julia_output/figures")
    fig = Figure(size=(1080, 450), fontsize=13)

    ax1 = Axis(fig[1, 1];
        title = "(A) True vs Parsimony Jumps across Reverse Zoonosis Rate ψ",
        xlabel = "Reverse zoonosis multiplier ψ (with ε=0.01, φ=4.0)",
        ylabel = "Expected jump count (n₁=n₂=100)"
    )
    lines!(ax1, df_psi.psi, df_psi.J_true_total; color=COL_BLACK, linewidth=2.6, label="True total (J*₁₂ + J*₂₁)")
    lines!(ax1, df_psi.psi, df_psi.J_true_12; color=COL_BLACK, linewidth=2.0, linestyle=:dot, label="True spillover J*₁₂")
    lines!(ax1, df_psi.psi, df_psi.J_true_21; color=COL_PURPLE, linewidth=2.0, linestyle=:dot, label="True reverse J*₂₁")
    lines!(ax1, df_psi.psi, df_psi.J_pars_1to2; color=COL_RES, linewidth=2.2, linestyle=:dash, label="Prior 1 (1→2 CS)")
    lines!(ax1, df_psi.psi, df_psi.J_pars_2to1; color=COL_HUM, linewidth=2.2, linestyle=:dashdot, label="Prior 2 (2→1 CS)")
    lines!(ax1, df_psi.psi, df_psi.J_pars_fitch; color=COL_REV, linewidth=2.2, label="Prior 3 (Fitch)")
    axislegend(ax1; position=:lt, framevisible=false, labelsize=10)

    ax2 = Axis(fig[1, 2];
        title = "(B) Parsimony Recovery Ratio Relative to Total Jumps J*",
        xlabel = "Reverse zoonosis multiplier ψ",
        ylabel = "Recovery ratio J_pars / (J*₁₂ + J*₂₁)"
    )
    lines!(ax2, df_psi.psi, df_psi.ratio_p1_total; color=COL_RES, linewidth=2.4, linestyle=:dash, label="Prior 1 / J*total")
    lines!(ax2, df_psi.psi, df_psi.ratio_p2_total; color=COL_HUM, linewidth=2.4, linestyle=:dashdot, label="Prior 2 / J*total")
    lines!(ax2, df_psi.psi, df_psi.ratio_fitch_total; color=COL_REV, linewidth=2.4, label="Prior 3 (Fitch) / J*total")
    axislegend(ax2; position=:rt, framevisible=false, labelsize=10)

    _save_both(fig, outdir, "fig9")
    return fig
end

"""
    generate_all_figures(; data_dir="julia_output/data", fig_dir="julia_output/figures", rng=StableRNG(2026))

Run the complete `Phylodynamics.jl` pipeline (ODEs, stochastic reaction trees, parsimony,
Mk2 ML, and structured coalescent profiles), save CSV tables to `data_dir`, and generate
all 9 figures (`fig1.pdf`–`fig9.pdf` + `.png`) in `fig_dir`.
"""
function generate_all_figures(;
    data_dir::AbstractString = "julia_output/data",
    fig_dir::AbstractString = "julia_output/figures",
    rng::AbstractRNG = StableRNG(2026)
)
    mkpath(data_dir)
    mkpath(fig_dir)

    # 1. Solve analytical coalescent ODEs for Models 1, 2, 3
    res1 = solve_coalescent_odes(model1_params())
    res2 = solve_coalescent_odes(model2_params())
    res3 = solve_coalescent_odes(model3_params())
    ode_summary = compute_canonical_summary()
    CSV.write(joinpath(data_dir, "ode_canonical_summary.csv"), ode_summary)

    # 2. Parameter sweeps (T_sample, n1/n2, R22, psi)
    df_time = sweep_sampling_time()
    sw_n = sweep_sampling_fraction()
    df_r22 = sweep_R22()
    df_psi = sweep_psi()
    CSV.write(joinpath(data_dir, "ode_sweep_sampling_time.csv"), df_time)
    CSV.write(joinpath(data_dir, "ode_sweep_n2.csv"), sw_n.df_n2)
    CSV.write(joinpath(data_dir, "ode_sweep_n1.csv"), sw_n.df_n1)
    CSV.write(joinpath(data_dir, "ode_sweep_R22.csv"), df_r22)
    CSV.write(joinpath(data_dir, "ode_sweep_psi.csv"), df_psi)

    # 3. Simulate trees from ModelingToolkit ReactionSystem (n1=100, n2=100)
    t1 = simulate_coalescent_tree(model1_params(); n1=100, n2=100, rng=rng)
    t2 = simulate_coalescent_tree(model2_params(); n1=100, n2=100, rng=rng)
    t3 = simulate_coalescent_tree(model3_params(); n1=100, n2=100, rng=rng)
    # Also save forward/ReMASTER stochastic trees
    sim1_stoch = simulate_remaster_tree(model1_params(); n1=100, n2=10, rng=rng)
    sim2_stoch = simulate_remaster_tree(model2_params(); n1=100, n2=10, rng=rng)
    sim3_stoch = simulate_remaster_tree(model3_params(); n1=100, n2=10, rng=rng)
    write_newick_file(sim1_stoch.binary_tree, joinpath(data_dir, "model1_remaster_stoch.nwk"); include_metadata=true)
    write_newick_file(sim2_stoch.binary_tree, joinpath(data_dir, "model2_remaster_stoch.nwk"); include_metadata=true)
    write_newick_file(sim3_stoch.binary_tree, joinpath(data_dir, "model3_remaster_stoch.nwk"); include_metadata=true)
    write_newick_file(t1, joinpath(data_dir, "model1_mtk_tree.nwk"); include_metadata=true)
    write_newick_file(t2, joinpath(data_dir, "model2_mtk_tree.nwk"); include_metadata=true)
    write_newick_file(t3, joinpath(data_dir, "model3_mtk_tree.nwk"); include_metadata=true)

    # 4. Fit parsimony and Mk2 ML on simulated ensembles
    sim_rows = NamedTuple[]
    for (mname, tr) in [("Model 1", t1), ("Model 2", t2), ("Model 3", t3)]
        pc = parsimony_counts(tr)
        mk = fit_mk2_ml(tr; n_maps=50, rng=rng)
        ch = cherry_assortativity(tr)
        sk = sackin_index(tr)
        push!(sim_rows, (
            model = mname,
            J_true_12 = Float64(tr.true_jumps_12),
            J_true_21 = Float64(tr.true_jumps_21),
            J_true_total = Float64(tr.true_jumps_12 + tr.true_jumps_21),
            J_pars_1to2 = Float64(pc.J_pars_1to2),
            J_pars_2to1 = Float64(pc.J_pars_2to1),
            J_pars_fitch = Float64(pc.J_pars_fitch),
            J_ml_unidir = mk.J_ml_unidir,
            J_ml_bidir_12 = mk.J_ml_bidir_12,
            J_ml_bidir_21 = mk.J_ml_bidir_21,
            J_ml_bidir_tot = mk.J_ml_bidir_tot,
            r_assort = ch.r,
            IS_norm = sk.IS_norm
        ))
    end
    df_sim_summary = DataFrame(sim_rows)
    CSV.write(joinpath(data_dir, "simulation_inference_summary.csv"), df_sim_summary)

    # 5. Structured coalescent profiles (PL2 engine)
    t1_prof = simulate_coalescent_tree(model1_params(); n1=100, n2=10, rng=StableRNG(101))
    t2_prof = simulate_coalescent_tree(model2_params(); n1=100, n2=10, rng=StableRNG(202))
    t3_prof = simulate_coalescent_tree(model3_params(); n1=100, n2=30, rng=StableRNG(42))
    prof = profile_structured_coalescent(t1_prof, t2_prof, t3_prof; likelihood=:PL2)
    CSV.write(joinpath(data_dir, "profile_epsilon_model1.csv"), prof.profile_eps)
    CSV.write(joinpath(data_dir, "profile_phi_model2.csv"), prof.profile_phi)
    CSV.write(joinpath(data_dir, "profile_psi_model3.csv"), prof.profile_psi)
    CSV.write(joinpath(data_dir, "surface_eps_phi_model2.csv"), prof.surface_2d)

    # 6. Render all 9 figures with CairoMakie
    plot_fig1_epidemic_and_coalescent(res1, res2, res3; outdir=fig_dir)
    plot_fig2_simulated_trees(t1, t2, t3; outdir=fig_dir)
    plot_fig3_sampling_time_bias(df_time; outdir=fig_dir)
    plot_fig4_sampling_fraction_bias(sw_n.df_n2, sw_n.df_n1; outdir=fig_dir)
    plot_fig5_R22_sweep(df_r22; outdir=fig_dir)
    plot_fig6_parsimony_vs_mk2(df_sim_summary, ode_summary; outdir=fig_dir)
    plot_fig7_structured_coalescent_profiles(prof; outdir=fig_dir)
    plot_fig8_topology_diagnostics(ode_summary, df_sim_summary; outdir=fig_dir)
    plot_fig9_reverse_zoonosis_sweep(df_psi; outdir=fig_dir)

    return (
        ode_summary = ode_summary,
        sim_summary = df_sim_summary,
        profiles = prof
    )
end
