using CartesianGrids
using Statistics: mean

Δx = 0.01
xlim = (-3,3)
ylim = (-3,3)
g = PhysicalGrid(xlim,ylim,Δx)

w = Nodes(Dual,size(g))
ψ = Nodes(Dual,w);
L = plan_laplacian(size(w),with_inverse=true)

xg,yg = coordinates(w,g);

c = 1.0
body = Plate(c,35)

# Find the minimum arc length
Δs = minimum(dlength(body))

# Move the airfoil
xc = 0.0; yc = 0.0
α = π/6
T = RigidTransform((xc,yc),-α)
T(body)

X = VectorData(midpoints(body)[1][1:end],midpoints(body)[2][1:end])
N = length(X.u)
f = ScalarData(X);
f̃ = ScalarData(X);
f₀ = ScalarData(X);
ψb = ScalarData(X);
ψ₀ = [0.0];

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
end
@testset "PotentialFlowSolution" begin
    ψ₀ = rand(1)
    δΓ_kvec = rand(1)
    sol = PotentialFlowSolution(ψ,f̃,ψ₀,δΓ_kvec)
    @test sol.ψ == ψ
    @test sol.f̃ == f̃
    @test sol.ψ₀ == ψ₀
    @test sol.δΓ_kvec == δΓ_kvec
    sol = PotentialFlowSolution(ψ,f̃,ψ₀)
    @test sol.ψ == ψ
    @test sol.f̃ == f̃
    @test sol.ψ₀ == ψ₀
    sol = PotentialFlowSolution(ψ,f̃)
    @test sol.ψ == ψ
    @test sol.f == f̃
end
