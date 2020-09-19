module GridPotentialFlow

using Reexport
using UnPack
using RecursiveArrayTools

@reexport using CartesianGrids
@reexport using RigidBodyTools
@reexport using ConstrainedSystems

include("bodyunitvectors.jl")
include("suctionparameter.jl")

include("solver/righthandside.jl")
include("solver/solution.jl")
include("solver/systems.jl")

include("vortex.jl")
include("vortexlist.jl")
include("vortexmodel.jl")

end
