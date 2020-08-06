module GridPotentialFlow

using Reexport
using UnPack
using RecursiveArrayTools

@reexport using CartesianGrids
@reexport using RigidBodyTools
@reexport using ConstrainedSystems

export BodyUnitVector, PotentialFlowSolution, PotentialFlowRHS, UnregularizedPotentialFlowSystem, RegularizedPotentialFlowSystem, Vortex, PotentialFlowSystem

include("vortex.jl")
include("bodyunitvectors.jl")
include("solver/righthandside.jl")
include("solver/solution.jl")
include("solver/systems.jl")


end
