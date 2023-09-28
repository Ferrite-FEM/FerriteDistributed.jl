using TimerOutputs

dto = TimerOutput()
reset_timer!(dto)

const liveserver = "liveserver" in ARGS

if liveserver
    using Revise
    @timeit dto "Revise.revise()" Revise.revise()
end

using Documenter
using FerriteDistributed

# Generate examples
include("generate.jl")

GENERATEDEXAMPLES = [joinpath("examples", f) for f in (
    "heat_equation_hypre.md",
    "heat_equation_pa.md",
)]

makedocs(
    sitename = "FerriteDistributed.jl",
    format = Documenter.HTML(),
    doctest = false,
    strict = false,
    draft = true,
    pages = Any[
        "Home" => "index.md",
        "Examples" => [GENERATEDEXAMPLES;],
        "Reference" => [
            "reference/interface.md",
            "reference/grid.md",
            "reference/utils.md",
        ]
    ]
)

if !liveserver
    @timeit dto "deploydocs" deploydocs(
        repo = "github.com/Ferrite-FEM/FerriteDistributed.jl.git",
        push_preview=true,
        versions = [
            "stable" => "v^",
            "v#.#",
            "dev" => "dev"
        ],
        forcepush = true,
    )
end

print_timer(dto)
