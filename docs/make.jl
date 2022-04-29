using CDLmotif
using Documenter

DocMeta.setdocmeta!(CDLmotif, :DocTestSetup, :(using CDLmotif); recursive=true)

makedocs(;
    modules=[CDLmotif],
    authors="Shane Kuei Hsien Chu (skchu@wustl.edu)",
    repo="https://github.com/kchu25/CDLmotif.jl/blob/{commit}{path}#{line}",
    sitename="CDLmotif.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://kchu25.github.io/CDLmotif.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/kchu25/CDLmotif.jl",
    devbranch="main",
)
