using Documenter
using GridPotentialFlow

makedocs(
    sitename = "GridPotentialFlow",
    format = Documenter.HTML(),
    modules = [GridPotentialFlow]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
