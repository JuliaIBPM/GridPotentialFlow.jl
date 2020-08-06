using GridPotentialFlow
using Test

@testset "System vectors" begin include("solver/systemvectors.jl") end
# @testset "Unregularized potential flow systems" begin include("unregularizedsystems.jl") end
@testset "Regularized potential flow systems" begin include("solver/regularizedsystems.jl") end
