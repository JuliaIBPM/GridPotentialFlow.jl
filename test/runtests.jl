using GridPotentialFlow
using Test

@testset "System vectors" begin include("systemvectors.jl") end
# @testset "Unregularized potential flow systems" begin include("unregularizedsystems.jl") end
@testset "Regularized potential flow systems" begin include("regularizedsystems.jl") end
