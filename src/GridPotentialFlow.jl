module GridPotentialFlow

using Reexport
using StructArrays
using DocStringExtensions
using RecipesBase

@reexport using CartesianGrids
@reexport using RigidBodyTools

include("suctionparameter.jl")
include("bodies.jl")

include("solver/righthandside.jl")
include("solver/solution.jl")
include("solver/systems.jl")

include("vortex.jl")
include("vortexmodel.jl")

end
