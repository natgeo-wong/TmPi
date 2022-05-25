using CreateTmPiDataset
using Documenter

DocMeta.setdocmeta!(CreateTmPiDataset, :DocTestSetup, :(using CreateTmPiDataset); recursive=true)

makedocs(;
    modules=[CreateTmPiDataset],
    authors="Nathanael Wong <natgeo.wong@outlook.com>",
    repo="https://github.com/natgeo-wong/CreateTmPiDataset.jl/blob/{commit}{path}#{line}",
    sitename="CreateTmPiDataset.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://natgeo-wong.github.io/CreateTmPiDataset.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/natgeo-wong/CreateTmPiDataset.jl",
    devbranch="main",
)
