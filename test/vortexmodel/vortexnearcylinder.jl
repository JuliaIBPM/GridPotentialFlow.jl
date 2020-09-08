
# grid
Δx = 0.05
xlim = (-2,2)
ylim = (-2,2)
g = PhysicalGrid(xlim,ylim,Δx);

# circle
Rc = 1
circle = Circle(Rc,Δx)

# vortices
v = Vortex(1.5,0.0,1.0);

@testset "Vortex near cylinder" begin
    # analytical solution
    Rv = sqrt(v.x^2+v.y^2)
    Vθ = -v.Γ/(2π*(Rv-Rc^2/Rv))
    Tv = 2π*Rv/abs(Vθ)

    # numerical solution
    Δt = 0.03
    T = 0:Δt:Tv
    X_hist = []
    model = VortexModel(g,bodies=circle,vortices=[v])

    for t in T
        Ẋ = computevortexvelocities(model)
        vortices = deepcopy(model.vortices.list)
        X = getvortexpositions(model)
        X = X + Ẋ*Δt
        setvortexpositions!(model,X)
        push!(X_hist,X)
    end

    @test isapprox(X_hist[end][1], 1.5; atol = 1e-1)
    @test isapprox(X_hist[end][2], 0.0; atol = 1e-1)
end
