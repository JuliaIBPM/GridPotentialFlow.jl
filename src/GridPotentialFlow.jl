module GridPotentialFlow

using Reexport
using UnPack
using StaticArrays
using StructArrays
using DocStringExtensions

@reexport using CartesianGrids
@reexport using RigidBodyTools

include("bodyunitvectors.jl")
include("modelparameters.jl")
include("bodies.jl")
include("utils.jl")

include("solver/righthandside.jl")
include("solver/solution.jl")
include("solver/systems.jl")

include("vortex.jl")
include("vortexlist.jl")
include("vortexmodel.jl")

end
