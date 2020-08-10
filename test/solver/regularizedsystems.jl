using LinearAlgebra

U∞ = 1.0

Δx = 0.0025
xlim = (-1,1)
ylim = (-1,1)
g = PhysicalGrid(xlim,ylim,Δx)

w = Nodes(Dual,size(g))
ψ = Nodes(Dual,w);
L = plan_laplacian(size(w),with_inverse=true)

xg,yg = coordinates(w,g);

c = 1.0
# body = Plate(c,Int(ceil(c/Δx)))
body = Plate(c,400)

# Find the minimum arc length
Δs = minimum(dlength(body))
# println("Ratio of arc spacing to cell size = ",Δs/Δx)

# Move the airfoil
xc = 0.0; yc = 0.0
α = π/3
T = RigidTransform((xc,yc),-α)
T(body)

# X = VectorData(midpoints(body.x)[1:end],midpoints(body.y)[1:end])
X = VectorData(body.x,body.y)
N = length(X.u)
f = ScalarData(X);
f̃ = ScalarData(X);
f₀ = ScalarData(X);
ψb = ScalarData(X);
ψ₀ = [0.0];

regop = Regularize(X,Δx,I0=origin(g),issymmetric=true)
Rmat,Emat = RegularizationMatrix(regop,f,w);
S = SaddleSystem(L,Emat,Rmat,SaddleVector(w,ψb))

# e1 = BodyUnitVector(N,1,RigidBodyTools.OpenBody);
# e2 = BodyUnitVector(N,N+1,RigidBodyTools.OpenBody);
e1 = BodyUnitVector(N,1,RigidBodyTools.OpenBody);
e2 = BodyUnitVector(N,N,RigidBodyTools.OpenBody);

ψb .= 1
w₀ = zero(w)
f₀ = constraint(S\SaddleVector(w₀,ψb));

Df₀ = Diagonal(f₀);
R̃mat = deepcopy(Rmat);
R̃mat.M .= R̃mat.M*Df₀;
S̃ = SaddleSystem(L,Emat,R̃mat,SaddleVector(w,ψb))

@testset "Internal functions" begin
    f̃ones = ScalarData(f₀)
    f̃ones .= 1.0
    f̃step = ScalarData(f₀)
    f̃step[1:Int(ceil(N/2))] .= -1.0
    f̃step[Int(ceil(N/2))+1:end] .= 1.0
    Nk = 2

    k_sheddingedges, activef̃limits_kvec = GridPotentialFlow._findsheddingedges(Nk, [e1,e2], f̃ones, [0.0,0.0], [0.0,0.0])
    @test k_sheddingedges == [1,2]
    @test activef̃limits_kvec == [0.0,0.0]

    k_sheddingedges, activef̃limits_kvec = GridPotentialFlow._findsheddingedges(Nk, [e1,e2], f̃ones, [-0.5,-0.5], [0.5,0.5])
    @test k_sheddingedges == [1,2]
    @test activef̃limits_kvec == [0.5,0.5]

    k_sheddingedges, activef̃limits_kvec = GridPotentialFlow._findsheddingedges(Nk, [e1,e2], -1*f̃ones, [-0.5,-0.5], [0.5,0.5])
    @test k_sheddingedges == [1,2]
    @test activef̃limits_kvec == [-0.5,-0.5]

    k_sheddingedges, activef̃limits_kvec = GridPotentialFlow._findsheddingedges(Nk, [e1,e2], f̃ones, [-2.0,-2.0], [2.0,2.0])
    @test isempty(k_sheddingedges)
    @test activef̃limits_kvec == [nothing,nothing]

    k_sheddingedges, activef̃limits_kvec = GridPotentialFlow._findsheddingedges(Nk, [e1,e2], f̃ones, [0.0,-2.0], [0.0,2.0])
    @test k_sheddingedges == [1]
    @test activef̃limits_kvec == [0.0,nothing]

    k_sheddingedges, activef̃limits_kvec = GridPotentialFlow._findsheddingedges(Nk, [e1,e2], f̃ones, [-2.0,0.0], [2.0,0.0])
    @test k_sheddingedges == [2]
    @test activef̃limits_kvec == [nothing,0.0]

    k_sheddingedges, activef̃limits_kvec = GridPotentialFlow._findsheddingedges(Nk, [e1,e2], f̃step, [0.0,0.0], [0.0,0.0])
    @test k_sheddingedges == [1,2]
    @test activef̃limits_kvec == [0.0,0.0]

    k_sheddingedges, activef̃limits_kvec = GridPotentialFlow._findsheddingedges(Nk, [e1,e2], f̃step, [-0.5,-0.5], [0.5,0.5])
    @test k_sheddingedges == [1,2]
    @test activef̃limits_kvec == [-0.5,0.5]

    #################################

    P_kvec = GridPotentialFlow._computesparsekuttaoperator.([e1,e2])
    δΓ_kvec = GridPotentialFlow._computepointvortexstrengths(Nk, Integer[], P_kvec, [f₀,f₀], [0.0,0.0], f₀, f₀, w)

    @test δΓ_kvec == [0.0, 0.0]

end

@testset "Steady (Kutta-Joukowski lift)" begin
    regularizedsys = PotentialFlowSystem(S̃,f₀,[e2]);
    regularizedsol = PotentialFlowSolution(ψ,f,[0.0],nothing)
    ψb .= -U∞*(X.v .- body.cent[2]);
    regularizedrhs = PotentialFlowRHS(w,ψb,[0.0],nothing,nothing)
    GridPotentialFlow.ldiv!(regularizedsol,regularizedsys,regularizedrhs)

    Γnumerical = sum(f₀.*regularizedsol.f̃)
    Γexact = -π*U∞*c*sin(α)

    @test isapprox(Γnumerical, Γexact; rtol = 1e-1)
end

@testset "Unsteady" begin
    Δt = 5e-2

    f̃ = constraint(S̃\SaddleVector(w,ψb))

    v1 = Vortex(body.x[1]-3Δt*c*cos(α+π/2),body.y[1]+3Δt*c*sin(α+π/2),1.0)
    d1 = Nodes(Dual,size(g));
    Rd1 =Regularize(VectorData([v1.x],[v1.y]),Δx,I0=origin(g),ddftype=CartesianGrids.Yang3)
    Rd1(d1,ScalarData([1.0]));

    v2 = Vortex(body.x[end]-3Δt*c*cos(α+π/2),body.y[end]+3Δt*c*sin(α+π/2),1.0)
    d2 = Nodes(Dual,size(g));
    Rd2 =Regularize(VectorData([v2.x],[v2.y]),Δx,I0=origin(g),ddftype=CartesianGrids.Yang3)
    Rd2(d2,ScalarData([1.0]));

    regularizedsys = PotentialFlowSystem(S̃,f₀,[e1,e2],[d1,d2]);
    regularizedsol = PotentialFlowSolution(ψ,f,[0.0],[0.0,0.0])
    ψb .= -U∞*(X.v .- body.cent[2]);
    regularizedrhs = PotentialFlowRHS(w,ψb,nothing,[0.0,0.0],[0.0,0.0])
    GridPotentialFlow.ldiv!(regularizedsol,regularizedsys,regularizedrhs)

    include("flatplatevalidation.jl")

    @test isapprox(regularizedsol.δΓ_kvec[2]/Δx^2, blobs[1].S; rtol = 1e-1) # edge indices are reversed

end
