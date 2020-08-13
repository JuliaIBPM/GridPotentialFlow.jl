module GridPotentialFlow

using Reexport
using UnPack
using RecursiveArrayTools

@reexport using CartesianGrids
@reexport using RigidBodyTools
@reexport using ConstrainedSystems

export BodyUnitVector, PotentialFlowSolution, PotentialFlowRHS, UnregularizedPotentialFlowSystem, RegularizedPotentialFlowSystem, PotentialFlowSystem, Vortex, updateposition!, VortexList, VortexModel, computeψ, computew, computevelocity, computeregularizationmatrix, getstrengths, getpositions, setvortexpositions!, getvortexpositions

include("bodyunitvectors.jl")
include("solver/righthandside.jl")
include("solver/solution.jl")
include("solver/systems.jl")

include("vortex.jl")
include("vortexlist.jl")
include("vortexmodel.jl")

end
