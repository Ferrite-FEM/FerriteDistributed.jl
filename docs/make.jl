using Documenter
using FerriteDistributed

makedocs(
    sitename = "FerriteDistributed.jl",
    format = Documenter.HTML(),
    doctest = false,
    strict = false,
    pages = Any[
        "Home" => "index.md",
        "API Reference" => [
            "reference/interface.md",
            "reference/grid.md",
            "reference/utils.md",
        ]
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
