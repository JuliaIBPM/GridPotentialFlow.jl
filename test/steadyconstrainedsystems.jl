using GridPotentialFlow
using Test

@testset "Steady constrained systems" begin
    Δx = 0.01
    xlim = (-1,1)
    ylim = (-2,2)
    g = PhysicalGrid(xlim,ylim,Δx)
    U∞ = 1.0
    params = Dict()
    params["freestream speed"] = U∞
    params["freestream angle"] = 0.0
    Δs = 1.4*cellsize(g)
    circle = Circle(0.5,Δs)
    plate = Plate(1.0,Δs)
    naca4 = NACA4(0.0,0.0,0.12,300)
    thickairfoil = SplinedBody(naca4.x,naca4.y,Δs)
    α = 10*π/180
    T1 = RigidTransform((0.0,1.0),0.0)
    T2 = RigidTransform((0.0,0.0),-α)
    T3 = RigidTransform((0.0,-1.0),-α)
    T1(circle)
    T2(plate)
    T3(thickairfoil)

    @testset "Circulation single open body" begin
        pfb = PotentialFlowBody(plate,Γ=1.0)
        prob = setup_problem(g,pfb,phys_params=params)
        sys = construct_system(prob)
        γ = vortexsheetstrength(sys)
        Γ = integrate(γ,sys.base_cache.ds)
        @test Γ ≈ 1.0
    end

    @testset "Circulation single closed body" begin
        pfb = PotentialFlowBody(circle,Γ=1.0)
        prob = setup_problem(g,pfb,phys_params=params)
        sys = construct_system(prob)
        γ = vortexsheetstrength(sys)
        Γ = integrate(γ,sys.base_cache.ds)
        @test Γ ≈ 1.0
    end

    @testset "Kutta condition single open body" begin
        pfb = PotentialFlowBody(plate,edges=[length(plate)+1])
        prob = setup_problem(g,pfb,phys_params=params)
        sys = construct_system(prob)
        γ = vortexsheetstrength(sys)
        Γ = integrate(γ,sys.base_cache.ds)
        @test isapprox(Γ,-π*U∞*sin(α),rtol=0.2)
        @test isapprox(-0.5*γ[end-1] + 1.5*γ[end],0.0,atol=1e-9)
    end

    @testset "Kutta condition single closed body" begin
        pfb = PotentialFlowBody(thickairfoil,edges=[length(thickairfoil)+1])
        prob = setup_problem(g,pfb,phys_params=params)
        sys = construct_system(prob)
        γ = vortexsheetstrength(sys)
        Γ = integrate(γ,sys.base_cache.ds)
        @test isapprox(Γ,-π*U∞*sin(α),rtol=0.2)
        γTE = 0.5*γ[1] + 0.5*γ[end]
        @test isapprox(γTE,0.0,atol=1e-9)
    end

    @testset "Circulation multibody" begin
        pfb1 = PotentialFlowBody(circle,Γ=-1.0)
        pfb2 = PotentialFlowBody(plate,Γ=1.0)
        pfb3 = PotentialFlowBody(thickairfoil,Γ=2.0)
        prob = setup_problem(g,BodyList([pfb1,pfb2,pfb3]),phys_params=params)
        sys = construct_system(prob)
        γ = vortexsheetstrength(sys)
        Γ1 = integrate(γ,sys.base_cache.ds,sys.base_cache.bl,1)
        Γ2 = integrate(γ,sys.base_cache.ds,sys.base_cache.bl,2)
        Γ3 = integrate(γ,sys.base_cache.ds,sys.base_cache.bl,3)
        @test Γ1 ≈ -1.0
        @test Γ2 ≈ 1.0
        @test Γ3 ≈ 2.0
    end

    @testset "Kutta condition multibody" begin
        pfb1 = PotentialFlowBody(plate,edges=[length(plate)+1])
        pfb2 = PotentialFlowBody(thickairfoil,edges=[length(thickairfoil)+1])
        prob = setup_problem(g,BodyList([pfb1,pfb2]),phys_params=params)
        sys = construct_system(prob)
        γ = vortexsheetstrength(sys)
        γTE1 = -0.5*γ[length(plate)-1] + 1.5*γ[length(plate)]
        γTE2 = 0.5*γ[length(plate)+1] + 0.5*γ[end]
        @test isapprox(γTE1,0.0,atol=1e-9)
        @test isapprox(γTE2,0.0,atol=1e-9)
    end

    @testset "Mixed multibody" begin
        pfb1 = PotentialFlowBody(circle,Γ=-1.0)
        pfb2 = PotentialFlowBody(plate,edges=[length(plate)+1])
        pfb3 = PotentialFlowBody(thickairfoil,edges=[length(thickairfoil)+1])
        prob = setup_problem(g,BodyList([pfb1,pfb2,pfb3]),phys_params=params)
        sys = construct_system(prob)
        γ = vortexsheetstrength(sys)
        Γ1 = integrate(γ,sys.base_cache.ds,sys.base_cache.bl,1)
        @test Γ1 ≈ -1.0
        γTE2 = -0.5*γ[length(circle)+length(plate)-1] + 1.5*γ[length(circle)+length(plate)]
        γTE3 = 0.5*γ[length(circle)+length(plate)+1] + 0.5*γ[end]
        @test isapprox(γTE2,0.0,atol=1e-9)
        @test isapprox(γTE3,0.0,atol=1e-9)
    end
end
