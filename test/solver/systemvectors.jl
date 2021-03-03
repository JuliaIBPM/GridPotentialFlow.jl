using CartesianGrids
using Statistics: mean


w = Nodes(Dual,(20,20))
ψ = Nodes(Dual,w);

w .= rand(20,20)
ψ .= rand(20,20)

c = 1.0
body = Plate(c,35)

X = VectorData(midpoints(body)[1][1:end],midpoints(body)[2][1:end])
N = length(X.u)
f = ScalarData(X);
f .= rand(N)
f₀ = ScalarData(X);
f₀ .= rand(N)
ψb = ScalarData(X);
ψb .= rand(N)
ψ₀ = rand(1,1);
Γw = randn(Float64)

@testset "BodyUnitVector" begin
    constantvector = ScalarData(N)
    constantvector .= 1.0
    step = 1.0
    linearlyincreasingvector = ScalarData(N)
    linearlyincreasingvector .= 1.0:step:N

    # Regular points
    e_1 = BodyUnitVector(N,1,RigidBodyTools.OpenBody,midpoints=false);
    e_end = BodyUnitVector(N,N,RigidBodyTools.OpenBody,midpoints=false);
    @test_throws AssertionError BodyUnitVector(N,N+1,RigidBodyTools.OpenBody,midpoints=false);
    @test_throws AssertionError BodyUnitVector(N,0,RigidBodyTools.OpenBody,midpoints=false);

    @test e_1'*constantvector == 1.0
    @test e_end'*constantvector == 1.0

    @test e_1'*linearlyincreasingvector == linearlyincreasingvector[1]
    @test e_end'*linearlyincreasingvector == linearlyincreasingvector[end]

    # Midpoints
    e_1 = BodyUnitVector(N,1,RigidBodyTools.OpenBody,midpoints=true);
    e_end = BodyUnitVector(N,N+1,RigidBodyTools.OpenBody,midpoints=true);
    @test_throws AssertionError BodyUnitVector(N,N+2,RigidBodyTools.OpenBody,midpoints=true);
    @test_throws AssertionError BodyUnitVector(N,0,RigidBodyTools.OpenBody,midpoints=true);

    @test e_1'*constantvector == 1.0
    @test e_end'*constantvector == 1.0

    @test e_1'*linearlyincreasingvector == linearlyincreasingvector[1]-step/2
    @test e_end'*linearlyincreasingvector == linearlyincreasingvector[end]+step/2
end
@testset "PotentialFlowRHS" begin
    rhs = PotentialFlowRHS(w,ψb)
    @test rhs.w == w
    @test rhs.ψb == ψb

    f̃lim_kvec = [SuctionParameter(randn(Float64))]
    rhs = PotentialFlowRHS(w,ψb,f̃lim_kvec)
    @test rhs.w == w
    @test rhs.ψb == ψb
    @test rhs.f̃lim_kvec == f̃lim_kvec

    f̃lim_kvec = [SuctionParameterRange(randn(Float64),randn(Float64))]
    @test f̃lim_kvec[1].min ≤ f̃lim_kvec[1].max
    rhs = PotentialFlowRHS(w,ψb,f̃lim_kvec)
    @test rhs.w == w
    @test rhs.ψb == ψb
    @test rhs.f̃lim_kvec == f̃lim_kvec

    f̃lim_kvec = [SuctionParameterRange(randn(Float64),randn(Float64))]
    @test f̃lim_kvec[1].min ≤ f̃lim_kvec[1].max
    rhs = PotentialFlowRHS(w,ψb,f̃lim_kvec,Γw)
    @test rhs.w == w
    @test rhs.ψb == ψb
    @test rhs.f̃lim_kvec == f̃lim_kvec
    @test rhs.Γw == Γw
end
@testset "PotentialFlowSolution" begin
    ψ₀ = rand(1)
    δΓ_kvec = rand(1)
    sol = PotentialFlowSolution(ψ,f,ψ₀,δΓ_kvec)
    @test sol.ψ == ψ
    @test sol.f == f
    @test sol.ψ₀ == ψ₀
    @test sol.δΓ_kvec == δΓ_kvec
    sol = PotentialFlowSolution(ψ,f,ψ₀)
    @test sol.ψ == ψ
    @test sol.f == f
    @test sol.ψ₀ == ψ₀
    sol = PotentialFlowSolution(ψ,f)
    @test sol.ψ == ψ
    @test sol.f == f
end
