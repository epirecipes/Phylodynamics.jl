using Documenter
using Phylodynamics

include("vignettes.jl")
tutorials = generate_vignette_pages(joinpath(@__DIR__, "..", "vignettes"),
    joinpath(@__DIR__, "src", "vignettes"))

makedocs(
    sitename = "Phylodynamics.jl",
    modules = [Phylodynamics],
    authors = "Simon D. W. Frost",
    doctest = true,
    checkdocs = :exports,
    format = Documenter.HTML(
        canonical = "https://epirecip.es/Phylodynamics.jl/",
        prettyurls = true,
    ),
    pages = [
        "Home" => "index.md",
        "Tutorials" => ["Overview" => "vignettes/index.md", tutorials.pages...],
        "Multi-deme coalescent equations" => "coalescent-moments.md",
        "API reference" => "api.md",
    ],
)

for alias in tutorials.aliases
    write(joinpath(@__DIR__, "build", alias),
        """<!doctype html><html lang="en"><head><meta charset="utf-8">
        <meta http-equiv="refresh" content="0;url=./">
        <title>Executable tutorial</title></head>
        <body><a href="./">Read the executable tutorial</a></body></html>""")
end
