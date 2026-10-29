using Test
using Phylodynamics
using Phylodynamics.Zoonosis
using Symbolics
using CSV
using DataFrames

function symbolic_zero(x)
    residual = Symbolics.simplify(Symbolics.simplify_fractions(Symbolics.expand(x)); expand=true)
    return isequal(residual, 0)
end

function sympy_expression(x)
    function encode(ex)
        ex isa Symbol && return string(ex)
        ex isa Integer && return string(ex)
        ex isa AbstractFloat && return string(ex)
        if ex isa Expr && ex.head === :call
            op = ex.args[1]
            args = encode.(ex.args[2:end])
            op === (+) && return "(" * join(args, "+") * ")"
            op === (*) && return "(" * join(args, "*") * ")"
            op === (/) && return "(" * join(args, "/") * ")"
            op === (^) && return "(" * join(args, "**") * ")"
            op === (-) && return length(args) == 1 ? "(-$(only(args)))" :
                "(" * join(args, "-") * ")"
        end
        error("Unsupported Symbolics expression in SymPy export: $ex")
    end
    return encode(Symbolics.toexpr(x))
end

@testset "Symbolics regeneration of the 37 two-host equations" begin
    z = [Symbolics.variable(Symbol("z$i")) for i in 1:37]
    Y1, Y2, X1, X2, beta, epsilon, phi, psi =
        @variables Y1 Y2 X1 X2 beta epsilon phi psi
    p = (; beta, epsilon, phi, psi)
    F = [beta * X1 * Y1 epsilon * beta * X2 * Y1;
         psi * beta * X1 * Y2 phi * beta * X2 * Y2]
    Y = [Y1, Y2]
    layout = Phylodynamics._coalescent_moment_layout(2,
        (:jumps, :coalescences, :branch_lengths, :cherries, :sackin, :fitch,
            :camin_sokal, :parsimony_lengths),
        [1, 2], nothing; max_states=1000, max_parsimony_classes=32)
    u = Num.(zeros(layout.nstates))
    u[layout.groups[:lineages]] .= z[1:2]
    u[layout.groups[:external]] .= z[3:6]
    u[layout.groups[:descendants]] .= z[10:11]
    for (name, indices) in ((:camin_sokal_12, [15 16; 17 18]),
                            (:camin_sokal_21, [19 20; 21 22]),
                            (:fitch, [23 24 25; 26 27 28]))
        H = reshape(view(u, layout.groups[name]), 2, :)
        program = layout.programs[name]
        for i in 1:2, a in 1:2
            H[i, program.leaf_classes[a]] = z[indices[i, a]]
        end
        if name === :fitch
            ambiguous = only(findall(message -> all(iszero, message), program.messages))
            H[:, ambiguous] .= z[indices[:, 3]]
        end
    end
    du = similar(u)
    Phylodynamics._coalescent_moment_rhs!(du, u, F, zeros(Int, 2, 2), Y, layout;
        descendant_mean=/)
    regenerated = Num.(zeros(37))
    regenerated[1:2] .= du[layout.groups[:lineages]]
    regenerated[3:6] .= du[layout.groups[:external]]
    regenerated[7:9] .= du[layout.groups[:cherries]]
    regenerated[10:11] .= du[layout.groups[:descendants]]
    regenerated[12] = only(du[layout.groups[:sackin]])
    F11_only = [F[1, 1] 0; 0 0]
    du11 = similar(u)
    Phylodynamics._coalescent_moment_rhs!(du11, u, F11_only, zeros(Int, 2, 2), Y, layout;
        descendant_mean=/)
    regenerated[13] = only(du11[layout.groups[:sackin]])
    regenerated[14] = regenerated[12] - regenerated[13]
    for (name, indices, score) in ((:camin_sokal_12, [15 16; 17 18], 31),
                                   (:camin_sokal_21, [19 20; 21 22], 32),
                                   (:fitch, [23 24 25; 26 27 28], 33))
        H = reshape(view(du, layout.groups[name]), 2, :)
        program = layout.programs[name]
        for i in 1:2, a in 1:2
            regenerated[indices[i, a]] = H[i, program.leaf_classes[a]]
        end
        if name === :fitch
            ambiguous = only(findall(message -> all(iszero, message), program.messages))
            regenerated[indices[:, 3]] .= H[:, ambiguous]
        end
        regenerated[score] = only(du[layout.groups[Symbol(name, :_score)]])
    end
    jumps = reshape(view(du, layout.groups[:jumps]), 2, 2)
    regenerated[29], regenerated[30] = jumps[1, 2], jumps[2, 1]
    lengths = du[layout.groups[:branch_lengths]]
    regenerated[34] = sum(lengths)
    regenerated[35], regenerated[37] = lengths
    regenerated[36] = du[layout.groups[:camin_sokal_12_lengths][
        layout.programs[:camin_sokal_12].leaf_classes[1]]]

    reference = similar(regenerated)
    Phylodynamics.Zoonosis._zoonosis_coalescent_rhs!(reference, z, p, [X1, Y1, X2, Y2];
        positive=identity, descendant_mean=/)
    rows = NamedTuple[]
    for index in 1:37
        substitutions = if index in 15:18
            Dict(z[1] => z[15] + z[16], z[2] => z[17] + z[18])
        elseif index in 19:22
            Dict(z[1] => z[19] + z[20], z[2] => z[21] + z[22])
        elseif index in 23:28
            Dict(z[1] => sum(z[23:25]), z[2] => sum(z[26:28]))
        else
            Dict{Num,Num}()
        end
        expression = substitute(regenerated[index], substitutions)
        expected = substitute(reference[index], substitutions)
        @testset "state $index" begin
            @test symbolic_zero(expression - expected)
        end
        push!(rows, (claim="legacy_$index",
            expression=sympy_expression(expression), expected=sympy_expression(expected)))
    end
    # Export an independently recognizable two-deme slice, including physical migration.
    H11, H12, H21, H22, F11, F12, F21, F22, G11, G12, G21, G22 =
        @variables H11 H12 H21 H22 F11 F12 F21 F22 G11 G12 G21 G22
    L11, L12, L21, L22, mass1, mass2 = @variables L11 L12 L21 L22 mass1 mass2
    F = [F11 F12; F21 F22]
    G = [G11 G12; G21 G22]
    u = Num.(zeros(layout.nstates))
    u[layout.groups[:lineages]] .= [H11 + H12, H21 + H22]
    u[layout.groups[:external]] .= [L11, L21, L12, L22]
    u[layout.groups[:descendants]] .= [mass1, mass2]
    H = reshape(view(u, layout.groups[:camin_sokal_12]), 2, :)
    H[:, layout.programs[:camin_sokal_12].leaf_classes] .= [H11 H12; H21 H22]
    Phylodynamics._coalescent_moment_rhs!(du, u, F, G, Y, layout; descendant_mean=/)
    dH = reshape(view(du, layout.groups[:camin_sokal_12]), 2, :)
    dL = reshape(view(du, layout.groups[:external]), 2, 2)
    independent = Dict(
        "A1" => du[layout.groups[:lineages][1]], "A2" => du[layout.groups[:lineages][2]],
        "H11" => dH[1, 1], "H12" => dH[1, 2], "H21" => dH[2, 1], "H22" => dH[2, 2],
        "J" => only(du[layout.groups[:camin_sokal_12_score]]),
        "X1" => substitute(du[layout.groups[:descendants][1]], Dict(mass1 => X1, mass2 => X2)),
        "X2" => substitute(du[layout.groups[:descendants][2]], Dict(mass1 => X1, mass2 => X2)),
        "L11" => dL[1, 1], "L12" => dL[1, 2], "L21" => dL[2, 1], "L22" => dL[2, 2],
        "C11" => du[layout.groups[:cherries][1]],
        "C12" => du[layout.groups[:cherries][2]],
        "C22" => du[layout.groups[:cherries][3]])
    for (claim, expression) in sort!(collect(independent); by=first)
        push!(rows, (claim=claim, expression=sympy_expression(expression),
            expected=sympy_expression(expression)))
    end
    if "--export" in ARGS
        output = joinpath(@__DIR__, "..", "verification", "generated-equations.csv")
        CSV.write(output, DataFrame(rows))
        println("Exported 37 regenerated equations and 16 event-kernel expressions to $output")
    end
end

@testset "Symbolic three-deme conservation" begin
    m = 3
    layout = Phylodynamics._coalescent_moment_layout(m,
        (:fitch, :sankoff, :coalescences, :sackin), [1, 2, 3],
        [0 1 2; 3 0 1; 2 3 0]; max_states=1000, max_parsimony_classes=32)
    Y = [Symbolics.variable(Symbol("Y$i")) for i in 1:m]
    F = [Symbolics.variable(Symbol("F$i$j")) for i in 1:m, j in 1:m]
    G = [Symbolics.variable(Symbol("G$i$j")) for i in 1:m, j in 1:m]
    for name in (:fitch, :sankoff)
        n = length(layout.programs[name].messages)
        u = Num.(zeros(layout.nstates))
        H = reshape(view(u, layout.groups[name]), m, n)
        H .= [Symbolics.variable(Symbol("H$i", :_, a)) for i in 1:m, a in 1:n]
        u[layout.groups[:lineages]] .= vec(sum(H; dims=2))
        u[layout.groups[:descendants]] .= [Symbolics.variable(Symbol("X$i")) for i in 1:m]
        du = similar(u)
        Phylodynamics._coalescent_moment_rhs!(du, u, F, G, Y, layout; descendant_mean=/)
        dH = reshape(view(du, layout.groups[name]), m, n)
        for i in 1:m
            @test symbolic_zero(sum(dH[i, :]) - du[layout.groups[:lineages][i]])
        end
        @test symbolic_zero(sum(du[layout.groups[:lineages]]) +
            sum(du[layout.groups[:coalescences]]))
        @test symbolic_zero(sum(du[layout.groups[:descendants]]))
    end
end
