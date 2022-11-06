using GridPotentialFlow
using Test

Δx = 0.01
xlim = (-1,1)
ylim = (-2,2)
g = PhysicalGrid(xlim,ylim,Δx)
Δs = 1.4*cellsize(g)
circle = Circle(0.5,Δs)
pfb = PotentialFlowBody(circle)

@testset "Uniform flow" begin
    U∞ = 1.0
    α = π/6
    params = Dict()
    params["freestream speed"] = U∞
    params["freestream angle"] = α

    prob = setup_problem(g,pfb,scaling=GridScaling,phys_params=params)
    sys = construct_system(prob)

    # Compute the scalar potential field and the corresponding velocity field
    ϕ = scalarpotential(sys)
    vel_ϕ = zeros_grid(sys);
    grad!(vel_ϕ,ϕ,sys);

    # Eldredge JCP 2022 eq 38
    v_df = zeros_grid(sys);
    regularize_normal!(v_df,sys.extra_cache.dϕtemp,sys)
    vel_ϕ .-= v_df;

    # Compute the streamfunction first from the velocity boundary conditions
    ψ_from_v_and_dv = streamfunction(sys);

    # Compute the streamfunction also from the jumps in scalar potential and normal velocity
    ψ_from_dϕ_and_dvn = zeros_gridcurl(sys);
    dv = zeros_surface(sys);
    dvn = ScalarData(dv);
    dψ = ScalarData(dv);
    nrm = normals(sys);
    prescribed_surface_jump!(dv,0.0,sys);
    dvn .= (nrm.u .* dv.u .+ nrm.v .* dv.v);
    solve!(ψ_from_dϕ_and_dvn,dψ,sys.extra_cache.dϕtemp,dvn,sys,0.0);
    ψ∞ = zeros_gridcurl(sys);
    vectorpotential_uniformvecfield!(ψ∞,U∞*cos(α),U∞*sin(α),sys.base_cache);
    ψ_from_dϕ_and_dvn .+= ψ∞;

    # Make sure both are approximately the same
    @test isapprox(ψ_from_v_and_dv,ψ_from_dϕ_and_dvn,rtol=1e-2)

    # Compute the corresponding velocity field
    vel_ψ = zeros_grid(sys);
    curl!(vel_ψ,ψ_from_v_and_dv,sys);

    # Make sure both velocity fields are approximately the same away from the circle
    @test isapprox(vel_ψ.u[2:end-1,g.N[2]÷6],vel_ϕ.u[2:end-1,g.N[2]÷6],rtol=1e-2)

    # Make sure both velocity fields provide an approximately zero normal velocity on the body
    vn_ϕ = ScalarData(dv);
    normal_interpolate!(vn_ϕ,vel_ϕ,sys)
    @test isapprox(vn_ϕ,ScalarData(dv),atol=1e-10)
    vn_ψ = ScalarData(dv);
    normal_interpolate!(vn_ψ,vel_ψ,sys)
    @test isapprox(vn_ψ,ScalarData(dv),atol=1e-1) # this one is less accurate
end

@testset "Body motion" begin
    U = 1.0
    α = π/6

    params = Dict()
    prob = setup_problem(g,
                         circle,
                         scaling=GridScaling,
                         phys_params=params,
                         motions=RigidBodyMotion(U*cos(α)+im*U*sin(α), 0.0))
    sys = construct_system(prob)

    # Compute the scalar potential field and the corresponding velocity field
    ϕ = scalarpotential(sys)
    vel_ϕ = zeros_grid(sys);
    grad!(vel_ϕ,ϕ,sys);

    # Eldredge JCP 2022 eq 38
    v_df = zeros_grid(sys);
    regularize_normal!(v_df,sys.extra_cache.dϕtemp,sys)
    vel_ϕ .-= v_df;

    # Compute the streamfunction first from the velocity boundary conditions
    ψ_from_v_and_dv = streamfunction(sys);

    # Compute the streamfunction also from the jumps in scalar potential and normal velocity
    ψ_from_dϕ_and_dvn = zeros_gridcurl(sys);
    dv = zeros_surface(sys);
    dvn = ScalarData(dv);
    dψ = ScalarData(dv);
    nrm = normals(sys);
    prescribed_surface_jump!(dv,0.0,sys);
    dvn .= (nrm.u .* dv.u .+ nrm.v .* dv.v);
    solve!(ψ_from_dϕ_and_dvn,dψ,sys.extra_cache.dϕtemp,dvn,sys,0.0);

    # Make sure both are approximately the same
    @test isapprox(ψ_from_v_and_dv,ψ_from_dϕ_and_dvn,rtol=1e-1)

    # Compute the corresponding velocity field
    vel_ψ = zeros_grid(sys);
    curl!(vel_ψ,ψ_from_v_and_dv,sys);

    # Make sure both velocity fields are approximately the same away from the circle
    @test isapprox(vel_ψ.u[2:end-1,g.N[2]÷6],vel_ϕ.u[2:end-1,g.N[2]÷6],rtol=1e-1)

    # Make sure both velocity fields provide a normal velocity on the body approximately equal to the prescribed normal velocity
    vb_prescribed = zeros_surface(sys)
    surface_velocity!(vb_prescribed, sys, 0.0)
    vbn_prescribed = (nrm.u .* vb_prescribed.u .+ nrm.v .* vb_prescribed.v);
    vn_ϕ = ScalarData(dv);
    normal_interpolate!(vn_ϕ,vel_ϕ,sys)
    @test isapprox(vn_ϕ,vbn_prescribed,atol=1e-10)
    vn_ψ = ScalarData(dv);
    normal_interpolate!(vn_ψ,vel_ψ,sys)
    @test isapprox(vn_ψ,vbn_prescribed,atol=1e-1) # this one is less accurate
end
