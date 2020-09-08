
# grid
Δx = 0.01
xlim = (-1,1)
ylim = (-1,1)
g = PhysicalGrid(xlim,ylim,Δx);

# vortices
v1 = Vortex(0.5,0.0,1.0);
v2 = Vortex(-0.5,0.0,1.0);

@testset "Vortex near cylinder" begin
    # analytical solution
    d = sqrt((v1.x-v2.x)^2+(v1.y-v2.y)^2)
    Vθ = v1.Γ/(2*π*d) # Analytical tangential velocity
    Tv = π*d/abs(Vθ) # Analytical period

    # numerical solution
    Δt = 0.01
    T = 0:Δt:Tv
    X_hist = []
    model = VortexModel(g,vortices=[v1,v2])

    for t in T
        Ẋ = computevortexvelocities(model)
        vortices = deepcopy(model.vortices.list)
        X = getvortexpositions(model)
        X = X + Ẋ*Δt
        setvortexpositions!(model,X)
        push!(X_hist,X)
    end

    @test isapprox(X_hist[end][1], 0.5; atol = 1e-1)
    @test isapprox(X_hist[end][2], -0.5; atol = 1e-1)
    @test isapprox(X_hist[end][3], 0.0; atol = 1e-1)
    @test isapprox(X_hist[end][4], 0.0; atol = 1e-1)
end
