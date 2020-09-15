using GridPotentialFlow
using Test

@testset "System vectors" begin
    include("solver/systemvectors.jl")
end
# @testset "Unregularized potential flow systems" begin
#     include("unregularizedsystems.jl")
# end
@testset "Regularized potential flow systems" begin
    include("solver/regularizedsystems.jl")
end
@testset "Vortex model possible flows" begin
    include("vortexmodel/possibleflows.jl")
end
@testset "Vortex dynamics" begin
    include("vortexmodel/vortexnearcylinder.jl")
    include("vortexmodel/corotatingvortices.jl")
end
@testset "Vortex shedding" begin
    include("vortexmodel/vortexshedding.jl")
end
@testset "Impulse" begin
    include("vortexmodel/cylinderimpulse.jl")
end
