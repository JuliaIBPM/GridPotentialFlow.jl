module GridPotentialFlow

using Reexport
using StructArrays
using DocStringExtensions
using RecipesBase

@reexport using CartesianGrids
@reexport using RigidBodyTools
@reexport using ImmersedLayers
using UnPack

export pressurejump!, pressure!, velocity!, surface_velocity!, sided_pressures

include("suctionparameter.jl")
include("bodies.jl")

include("solver/problem.jl")
include("solver/righthandside.jl")
include("solver/solution.jl")
include("solver/systems.jl")

include("vortex.jl")
include("vortexmodel.jl")
include("pressure.jl")


end
