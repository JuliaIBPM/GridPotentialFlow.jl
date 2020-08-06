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
    e_1 = BodyUnitVector(N,1,RigidBodyTools.OpenBody);
    e_mid = BodyUnitVector(N,Int(ceil((N+1)/2)),RigidBodyTools.OpenBody);
    e_end = BodyUnitVector(N,N+1,RigidBodyTools.OpenBody);
    @test_throws AssertionError BodyUnitVector(N,N+2,RigidBodyTools.OpenBody);
    @test_throws AssertionError BodyUnitVector(N,0,RigidBodyTools.OpenBody);

    constantvector = ScalarData(N)
    constantvector .= 1.0
    step = 1.0
    linearlyincreasingvector = ScalarData(N)
    linearlyincreasingvector .= 1.0:step:N

    @test e_1'*constantvector == 1.0
    @test e_mid'*constantvector == 1.0
    @test e_end'*constantvector == 1.0

    @test e_1'*linearlyincreasingvector == linearlyincreasingvector[1]-step/2
    @test e_mid'*linearlyincreasingvector == mean(linearlyincreasingvector)
    @test e_end'*linearlyincreasingvector == linearlyincreasingvector[end]+step/2
end
@testset "PotentialFlowRHS" begin
    f̃limit_kvec = rand(1)
    steadyrhs = PotentialFlowRHS(w,ψb,f̃limit_kvec,nothing,nothing)
    @test steadyrhs.w == w
    @test steadyrhs.ψb == ψb
    @test steadyrhs.f̃limit_kvec == f̃limit_kvec
    @test steadyrhs.f̃min_kvec == nothing
    @test steadyrhs.f̃max_kvec == nothing
    unsteadyrhs = PotentialFlowRHS(w,ψb,nothing,[0.0],[0.0])
    @test unsteadyrhs.w == w
    @test unsteadyrhs.ψb == ψb
    @test unsteadyrhs.f̃limit_kvec == nothing
    @test unsteadyrhs.f̃min_kvec == [0.0]
    @test unsteadyrhs.f̃max_kvec == [0.0]
    rhs = PotentialFlowRHS(w,ψb)
    @test rhs.w == w
    @test rhs.ψb == ψb
end
@testset "PotentialFlowSolution" begin
    ψ₀ = rand(1)
    δΓ_kvec = rand(1)
    sol = PotentialFlowSolution(ψ,f̃,ψ₀,δΓ_kvec)
    @test sol.ψ == ψ
    @test sol.f̃ == f̃
    @test sol.ψ₀ == ψ₀
    @test sol.δΓ_kvec == δΓ_kvec
    sol = PotentialFlowSolution(ψ,f̃)
    @test sol.ψ == ψ
    @test sol.f̃ == f̃
end
