module GridPotentialFlow

using Reexport
using UnPack
using RecursiveArrayTools

@reexport using CartesianGrids
@reexport using RigidBodyTools
@reexport using ConstrainedSystems

export BodyUnitVector, PotentialFlowSolution, PotentialFlowRHS, UnregularizedPotentialFlowSystem, RegularizedPotentialFlowSystem, Vortex, PotentialFlowSystem

include("bodyunitvectors.jl")
include("righthandside.jl")
include("solution.jl")
include("vortices.jl")
include("systems.jl")


end
