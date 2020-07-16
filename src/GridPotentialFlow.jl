module GridPotentialFlow

using Reexport
using UnPack
using RecursiveArrayTools

@reexport using CartesianGrids
@reexport using RigidBodyTools
@reexport using ConstrainedSystems

export BodyUnitVector, SolutionVector, RightHandSideVector, UnregularizedPotentialFlowSystem, RegularizedPotentialFlowSystem, Vortex, PotentialFlowSystem

include("bodyunitvectors.jl")
include("systemvectors.jl")
include("vortices.jl")
include("systems.jl")


end
