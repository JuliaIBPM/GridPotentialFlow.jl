using GridPotentialFlow
using Plots
Δx = 0.02
Lx = 4.0
xlim = (-Lx/2,Lx/2)
ylim = (-Lx/2,Lx/2)
g = PhysicalGrid(xlim,ylim,Δx);

Rc = Lx/4
Δs = 2Δx
body = PotentialFlowBody(Circle(Rc,Δs))

Rv = 3/2*Rc
Γv = 1.0
v = Vortex(Rv,0.0,Γv);

setΓ(body,-Γv)

model = VortexModel(g,vortices=[v],bodies=[body]);
sol = solve(model);

r1 = sqrt.((body.points.x.-Rv).^2 .+ (body.points.y).^2) # distance from vortex to points on circle
r2 = sqrt.((body.points.x.-Rc^2/Rv).^2 .+ (body.points.y).^2) # distance from image vortex to points on circle
θ1 = π.+atan.((body.points.y)./(body.points.x.-Rv)) # angle from vortex to points on circle
θ2 = atan.((body.points.y),(body.points.x.-Rc^2/Rv)) # angle from image vortex to points on circle
v1x = -v.Γ./(2π*r1).*cos.(θ1.-π/2) # x velocity on circle induced by vortex
v1y = -v.Γ./(2π*r1).*sin.(θ1.-π/2) # y velocity on circle induced by vortex
v2x = v.Γ./(2π*r2).*cos.(θ2.-π/2) # x velocity on circle induced by image vortex
v2y = v.Γ./(2π*r2).*sin.(θ2.-π/2) # y velocity on circle induced by image vortex
V = sqrt.((v1x.+v2x).^2+(v1y.+v2y).^2) # velocity magnitude on circle
γ = -V; # bound vortex sheet strength on circle (velocity on circle is clockwise if positive vortex is to the right of it)

plot(sol.ψ,g)
plot!(body,fillcolor=:black,fillrange=0,fillalpha=0.25,linecolor=:black,linewidth=2)
scatter!([v.x],[v.y],color=:black,markersize=2,xlabel="x",ylabel="y")

plot(sol.f./Δs,label="f/ds",xlabel="body point index")
plot!(γ,label="gamma")

Vv = Γv/(2π*(Rv-Rc^2/Rv))
T = 2π*Rv/Vv

import OrdinaryDiffEq
function rhs(X,model,t)
    setvortexpositions!(model,X)
    Ẋ = vortexvelocities!(model)
    return Ẋ
end
X = getvortexpositions(model)
prob = OrdinaryDiffEq.ODEProblem(rhs,X,(0.0,T),model);
sol = OrdinaryDiffEq.solve(prob,dt=0.1,OrdinaryDiffEq.RK4(),dense=false,adaptive=false);
plot(body,fillcolor=:black,fillrange=0,fillalpha=0.25,linecolor=:black,linewidth=2)
plot!(map(s->s.u[1],sol.u),map(s->s.v[1],sol.u))
scatter!([v.x],[v.y],color=:black,markersize=2,xlabel="x",ylabel="y")


@testset "Vortex near cylinder" begin
    import OrdinaryDiffEq
#jl

    Rv = sqrt(v.x^2+v.y^2)
    Vθ = -v.Γ/(2π*(Rv-Rc^2/Rv))
    Tv = 2π*Rv/abs(Vθ)
#jl

    tspan = (0,Tv)
    model = VortexModel(g,bodies=[body],vortices=[v])
#jl
    function rhs(X,model,t)
        setvortexpositions!(model,X)
        Ẋ = vortexvelocities!(model)
        return Ẋ
    end
#jl
    X = getvortexpositions(model)
    prob = OrdinaryDiffEq.ODEProblem(rhs,X,tspan,model);
    sol = OrdinaryDiffEq.solve(prob,dt=0.1,OrdinaryDiffEq.RK4(),dense=false,adaptive=false);
#jl
    @test isapprox(sol.u[end].u[1], 1.5; atol = 1e-2)
    @test isapprox(sol.u[end].v[1], 0.0; atol = 1e-1)
end

body = PotentialFlowBody(Circle(Rc,Δs), U=(1.0,0.0))
model = VortexModel(g,bodies=[body]);
sol = solve(model);

θ = atan.(body.points.y,body.points.x)
γ = 2*sin.(θ);

plot(sol.ψ,g)
plot!(body,fillcolor=:black,fillrange=0,fillalpha=0.25,linecolor=:black,linewidth=2,xlabel="x",ylabel="y")

plot(sol.f./Δs,label="f/ds",xlabel="body point index")
plot!(γ,label="gamma")

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
