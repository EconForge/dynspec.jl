using Documenter
using Dynspec

makedocs(
    sitename = "Dynspec.jl",
    authors = "Pablo Winant and EconForge contributors",
    repo = Documenter.Remotes.GitHub("EconForge", "dynspec.jl"),
    checkdocs = :exports,
    warnonly = [:missing_docs],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        sidebar_sitename = true
    ),
    modules = [Dynspec, Dynspec.Dyno],
    pages = [
        "Home" => "index.md",
        "Function & Option Guide" => "guide.md",
        "API Reference" => "api.md",
    ]
)

deploydocs(
    repo = "github.com/EconForge/dynspec.jl.git",
    devbranch = "main",
)

