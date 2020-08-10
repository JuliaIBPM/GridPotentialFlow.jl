# Streamfunction on the body
Ub = 0.0
Vb = 0.0
Ωb = 0.0
U∞ = 1.0
V∞ = 0.0

# Body parameters
a = 0.5
b = 0.2

# Build high-res ellipse
n = 100;
body = Ellipse(a,b,n)

# Move the ellipse
xc = 0.0; yc = 0.0
T = RigidTransform((xc,yc),0.0)
T(body)

Xhires = VectorData(body.x,body.y)
N = length(Xhires.u);

lx = 2.0
xlim = (-lx/2,lx/2)
ylim = (-lx/2,lx/2)

error_hist = []
nx_hist = [2^k for k in 7:8]

for nx in nx_hist

    Δx = lx/nx
    g = PhysicalGrid(xlim,ylim,Δx)

    w = Nodes(Dual,size(g))
    xg,yg = coordinates(w,g);

    ψ = Nodes(Dual,size(g));

    L = plan_laplacian(size(ψ),with_inverse=true)
    L⁻¹(ψ::T) where {T} = L\ψ;

    bodySpline = RigidBodyTools.SplinedBody(hcat(Xhires.u[:],Xhires.v[:]),2.0Δx)

    Δs = minimum(Bodies.dlength(bodySpline))
    X = VectorData(bodySpline.x,bodySpline.y)
    N = length(X.u)

    f = ScalarData(N);
    One = ScalarData(N)
    One .= 1.0

    ψb = ScalarData(N)
    ψb .= Ub*X.v - Vb*X.u - 1/2*Ωb*(X.u .* X.u + X.v .* X.v) - (U∞*X.v - V∞*X.u);

    regop = Regularize(X,Δx,I0=origin(g),issymmetric=true,weights=Δs,ddftype=Fields.Yang3)
    Rmat,_ = RegularizationMatrix(regop,f,w);
    Emat = InterpolationMatrix(regop,w,f);
    S = SaddleSystem((ψ,f),(L⁻¹,Rmat,Emat),issymmetric=true)

    unregularizedsys = PotentialFlowSystem(S)
    unregularizedsol = PotentialFlowSolution(ψ,f,[0.0],[0.0])
    unregularizedrhs = PotentialFlowRHS(w,ψb)
    GridPotentialFlow.ldiv!(unregularizedsol,unregularizedsys,unregularizedRHS)

    γ = -U∞*(a+b)*a*X.v./sqrt.(b^4*X.u.^2+a^4*X.v.^2);
    error = norm(unregularizedsol.f./Δs.-γ)/N
    push!(error_hist,error)

end

println(error_hist)

@test log(error_hist[end-1]/error_hist[end])/log(2) > 1
