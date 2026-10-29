using Documenter
using Phylodynamics

vignettes = joinpath(@__DIR__, "..", "vignettes")
for (directory, _, files) in walkdir(vignettes)
    for file in files
        endswith(file, ".qmd") || continue
        html = replace(file, r"\.qmd$" => ".html")
        source = joinpath(directory, html)
        isfile(source) || error("Render the vignette before building documentation: $source")
        destination = joinpath(@__DIR__, "src", "vignettes",
            relpath(directory, vignettes), html)
        mkpath(dirname(destination))
        cp(source, destination; force = true)
    end
end

makedocs(
    sitename = "Phylodynamics.jl",
    modules = [Phylodynamics],
    authors = "Simon D. W. Frost",
    doctest = false,
    checkdocs = :none,
    format = Documenter.HTML(
        canonical = "https://epirecip.es/Phylodynamics.jl/",
        prettyurls = true,
    ),
    pages = [
        "Home" => "index.md",
        "Vignettes" => "vignettes/index.md",
        "API reference" => "api.md",
    ],
)
