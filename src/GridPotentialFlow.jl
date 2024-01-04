module GridPotentialFlow

using Reexport
using StructArrays
using DocStringExtensions
using RecipesBase
using UnPack

@reexport using CartesianGrids
@reexport using RigidBodyTools
@reexport using ImmersedLayers

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
