using TmPi
using Documenter

DocMeta.setdocmeta!(TmPi, :DocTestSetup, :(using TmPi); recursive=true)

makedocs(;
    modules=[TmPi],
    authors="Nathanael Wong <natgeo.wong@outlook.com>",
    repo="https://github.com/natgeo-wong/TmPi.jl/blob/{commit}{path}#{line}",
    sitename="TmPi.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://natgeo-wong.github.io/TmPi.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/natgeo-wong/TmPi.jl",
    devbranch="main",
)
