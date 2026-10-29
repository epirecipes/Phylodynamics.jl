module DocumentationConversionTests
using Test
include(joinpath(@__DIR__, "..", "docs", "vignettes.jl"))

@testset "Executable tutorial conversion" begin
    source = raw"""---
title: "A test tutorial"
---

An inline expression $x + y$ and a display expression $$x^2$$.

```{julia}
x = 1
println("$x")
```

```{julia}
@assert x == 1
x + 1
```

```sh
echo "$UNCHANGED"
```
"""
    page = documenter_vignette(source, "shared")
    @test startswith(page, "# A test tutorial\n")
    @test length(findall("```@example shared", page)) == 2
    @test occursin("``x + y``", page)
    @test occursin("```math\nx^2\n```", page)
    @test occursin(raw"""println("$x")""", page)
    @test occursin(raw"$UNCHANGED", page)
    @test_throws ErrorException documenter_vignette("No frontmatter", "bad")
    @test_throws ErrorException documenter_vignette(
        "---\ntitle: \"Missing fence\"\n---\n```{julia}\nx = 1\n", "bad")
    root = joinpath(@__DIR__, "..", "vignettes")
    mktempdir() do directory
        generated = generate_vignette_pages(root, directory)
        expected = 0
        for folder in filter(path -> isdir(joinpath(root, path)), readdir(root))
            for file in filter(path -> endswith(path, ".qmd"), readdir(joinpath(root, folder)))
                original = read(joinpath(root, folder, file), String)
                converted = read(joinpath(directory, folder, "index.md"), String)
                @test length(findall("```@example ", converted)) ==
                    length(findall("```{julia}", original))
                @test !occursin("```{julia}", converted)
                @test "vignettes/$folder/index.md" in last.(generated.pages)
                expected += 1
            end
        end
        @test length(generated.pages) == length(generated.aliases) == expected
    end
end
end
