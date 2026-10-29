function _documenter_math(text)
    text = replace(text, r"\$\$(.*?)\$\$"s =>
        value -> "\n```math\n$(strip(chop(value; head=2, tail=2)))\n```\n")
    return replace(text, r"(?<!\\)\$(?!\$)([^$\n]+?)(?<!\\)\$(?!\$)" =>
        value -> "``$(chop(value; head=1, tail=1))``")
end

function documenter_vignette(source, name)
    frontmatter = match(r"\A---\r?\n(.*?)\r?\n---\r?\n"s, source)
    frontmatter === nothing && error("Vignette $name is missing YAML frontmatter")
    title = match(r"(?m)^title: \"(.+)\"\r?$", frontmatter.captures[1])
    title === nothing && error("Vignette $name requires a quoted title")
    body = chop(source; head=length(frontmatter.match), tail=0)
    output = IOBuffer()
    text = IOBuffer()
    println(output, "# ", title.captures[1], "\n")
    println(output, "Every Julia block on this page is executed during the documentation build.\n")
    fence = :none
    examples = 0
    function flush_text!()
        print(output, _documenter_math(String(take!(text))))
    end
    for line in eachline(IOBuffer(body); keep=true)
        if fence === :none && occursin(r"^```\{julia\}\s*$", line)
            flush_text!()
            println(output, "```@example ", name)
            fence = :example
            examples += 1
        elseif fence === :none && startswith(line, "```")
            flush_text!()
            print(output, line)
            fence = :other
        elseif fence !== :none && occursin(r"^```\s*$", line)
            print(output, line)
            fence = :none
        elseif fence !== :none
            print(output, line)
        else
            print(text, line)
        end
    end
    fence === :none || error("Vignette $name has an unclosed code fence")
    flush_text!()
    examples > 0 || error("Vignette $name has no executable Julia blocks")
    return String(take!(output))
end

function generate_vignette_pages(vignettes, destination)
    pages = Pair{String,String}[]
    aliases = String[]
    for directory in sort(readdir(vignettes))
        source_dir = joinpath(vignettes, directory)
        isdir(source_dir) || continue
        sources = filter(file -> endswith(file, ".qmd"), readdir(source_dir))
        isempty(sources) && continue
        length(sources) == 1 || error("Expected one Quarto source in $source_dir")
        filename = only(sources)
        name = "vignette_" * replace(directory, r"[^A-Za-z0-9_]" => "_")
        markdown = documenter_vignette(read(joinpath(source_dir, filename), String), name)
        target = joinpath(destination, directory)
        mkpath(target)
        write(joinpath(target, "index.md"), markdown)
        old_html = replace(filename, r"\.qmd$" => ".html")
        # Remove the old copied artifact: only freshly executed output is published.
        rm(joinpath(target, old_html); force=true)
        push!(aliases, joinpath("vignettes", directory, old_html))
        label = titlecase(replace(directory, r"^\d+-" => "", "-" => " "))
        push!(pages, label => joinpath("vignettes", directory, "index.md"))
    end
    isempty(pages) && error("No executable vignettes found in $vignettes")
    return (; pages, aliases)
end
