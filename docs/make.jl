using Documenter
using GridPotentialFlow

makedocs(
    sitename = "GridPotentialFlow",
    format = Documenter.HTML(),
    modules = [GridPotentialFlow]
)
makedocs(
    sitename = "GridPotentialFlow.jl",
    doctest = true,
    clean = true,
    pages = [
        "Home" => "index.md",
        "Manual" => ["manual/1.-Basic-potential-flow-problem.md",
                     "manual/2.-Potential-flow-with-an-impenetrable-surface.md",
                     "manual/3.-Non-uniqueness-and-discrete-circulation.md",
                     "manual/4.-The-Kutta-condition.md",
                     "manual/5.-Generalized-edge-conditions.md",
                     "manual/6.-Force-and-the-added-mass.md"
                     ],
        "Library" => ["public.md",
                      "private.md"]
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        mathengine = MathJax2(Dict(
            :TeX => Dict(
                :equationNumbers => Dict(:autoNumber => "AMS"),
                :Macros => Dict()
            )
        ))
    )
)

deploydocs(
    repo = "github.com/JuliaIBPM/GridPotentialFlow.jl.git",
    target = "build"
)
