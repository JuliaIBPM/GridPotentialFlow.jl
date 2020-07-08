module GridPotentialFlow

using Reexport
using UnPack

@reexport using CartesianGrids
@reexport using RigidBodyTools
@reexport using ConstrainedSystems

export PotentialFlowSystem, Vortex

include("systems.jl")
include("vortices.jl")

end
