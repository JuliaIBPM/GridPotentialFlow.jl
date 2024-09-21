using GridPotentialFlow
Δx = 0.01
Lx = 2.0
xlim = (-Lx/2,Lx/2)
ylim = (-Lx/2,Lx/2)
g = PhysicalGrid(xlim,ylim,Δx,optimize=false);

v = Vortex(0.0,0.0,1.0);

model = VortexModel(g,vortices=[v]);

s = streamfunction(model);
using Plots
plot(s,g,xlabel="x",ylabel="y")

q = curl(s);
plot(q,g,xlabel="x",ylabel="y")

function ψ_vortex!(ψ::Nodes{Dual,nx,ny},vortex::Vortex,g::PhysicalGrid) where {nx,ny}
    x,y = coordinates(ψ,g)
    for i in 2:nx-1, j in 2:ny-1
        r = sqrt((x[i]-vortex.x)^2+(y[j]-vortex.y)^2)
        ψ[i,j] = ψ_vortex(r,vortex.Γ)
    end
end

function ψ_vortex(r::Real,Γ::Real)
    return -Γ/(2π)*log(r)
end

nv = 4;
vl = Vortex.(-Lx/4 .+ 0.4*Lx*(rand(nv).-0.5),-Lx/4 .+ 0.4*Lx*(rand(nv).-0.5),0.5*rand(nv).+0.5);

grids = [PhysicalGrid(xlim,ylim,Lx/(nx-2),optimize=false) for nx in [2^p for p in 5:9]];

using LinearAlgebra: norm

errors = []
gridspacings = []
for g in grids
    # LGF
    local model = VortexModel(g,vortices=vl)
    s_lgf = streamfunction(model);

    # Exact solution
    s_exact = Nodes(Dual,size(g))
    s_temp = Nodes(Dual,size(g))
    for v in vl
        ψ_vortex!(s_temp,v,g)
        s_exact += s_temp
    end

    # Bring to same reference level. A constant value can be added/subtracted to any potential flow solution
    s_lgf .-= s_lgf[g.I0[1],2];
    s_exact .-= s_exact[g.I0[1],2];

    error = s_lgf-s_exact
    idx = g.I0[1] + 1 : g.N[1] - 1 # Only look at top right corner
    push!(errors,norm(error[idx,idx])/norm(s_lgf[idx,idx]))
    push!(gridspacings,g.Δx)
end

firstorder=(errors[end]/gridspacings[end].*gridspacings).^1
secondorder=(sqrt(errors[end])/gridspacings[end].*gridspacings).^2
p=plot(gridspacings,errors,xaxis=:log,yaxis=:log,marker=:circle,lab="error",xlabel="dx",legend=:bottomright,title="Error")
plot!(gridspacings,firstorder,lab="1st order",linestyle=:dot,linecolor=:black)
plot!(gridspacings,secondorder,lab="2nd order",linestyle=:dash,linecolor=:black)

d = Lx/2
Γ = 1
v1 = Vortex(d/2,0.0,Γ);
v2 = Vortex(-d/2,0.0,Γ);

Vθ = Γ/(2*π*d); # Analytical tangential velocity
T = π*d/Vθ # Analytical period

model = VortexModel(g,vortices=[v1,v2]);

tspan = (0.0,T);

import OrdinaryDiffEq

#= To construct the `ODEProblem` from this package, we will have to provide a right-hand side function that computes the velocities $Ẋ$ of the vortices, which is the local flow velocity at their positions. These velocities are obtained using `vortexvelocities`, which regularizes the vorticity to the grid, solves the potential flow system, and interpolates the velocities from the grid to the vortex locations as

$\left(U_{q}, V_{q}\right)=\sum_{i, j} \mathsf{v}_{i j} d\left(\frac{\mathsf{x}_{i}-X_{q}}{\Delta x}\right)\left(\frac{\mathsf{y}_{j}-Y_{q}}{\Delta x}\right),$

where $v$ is the velocity field on the nodes.
=#

function rhs(X,model,t)
    setvortexpositions!(model,X)
    Ẋ = vortexvelocities!(model)
    return Ẋ
end

X = getvortexpositions(model)
prob = OrdinaryDiffEq.ODEProblem(rhs,X,tspan,model);
sol = OrdinaryDiffEq.solve(prob,dt=0.1,OrdinaryDiffEq.RK4(),dense=false,adaptive=false);

plot(sol.t,map(s->s.u[1],sol.u),xlabel='t',label="x₁")
plot!(sol.t,map(s->s.v[1],sol.u),label="y₁")

@testset "Corotating point vortices" begin
    @test isapprox(sol.u[end].u[1], 0.5; atol = 1e-2)
    @test isapprox(sol.u[end].u[2], -0.5; atol = 1e-2)
    @test isapprox(sol.u[end].v[1], 0.0; atol = 1e-2)
    @test isapprox(sol.u[end].v[2], 0.0; atol = 1e-2)
end

xlim = (-2,2);
ylim = (-2,2);
Δx = 0.05;
g = PhysicalGrid(xlim,ylim,Δx,optimize=false);

function vortexpatch!(vort,xc,yc,Γ,radius,nring)
    Δr = radius/(nring-1/2)
    dΓ = Γ/(1+8*nring*(nring-1)/2)
    push!(vort,Vortex(xc,yc,dΓ))
    for ir in 1:nring-1
        nθ = 8*ir
        for j = 0:nθ-1
            push!(vort,Vortex(xc + ir*Δr*cos(2π*j/nθ),yc + ir*Δr*sin(2π*j/nθ),dΓ))
        end
    end
    return vort
end

vortexpatch(xc,yc,Γ,radius,nring) = vortexpatch!(Vortex[],xc,yc,Γ,radius,nring)

d0 = d
anims = []
for r0 in [0.2*d0, 0.4*d0]
    vortices = vcat(vortexpatch(0.0,0.0+d0/2,Γ,r0,10),vortexpatch(0.0,0.0-d0/2,Γ,r0,10));
    vortexcolors = vcat(fill(:blue,length(vortices)÷2),fill(:red,length(vortices)÷2));
    local model = VortexModel(g,vortices=vortices);
    local X = getvortexpositions(model)
    local prob = OrdinaryDiffEq.ODEProblem(rhs,X,tspan,model);
    local sol = OrdinaryDiffEq.solve(prob,dt=0.1,OrdinaryDiffEq.RK4(),dense=false,adaptive=false);
    anim = @animate for i=1:length(sol.t)-1
        plot(xlims=xlim,ylims=ylim,ratio=:equal,legend=:none,title="r0/d0 = $(r0/d0)")
        scatter!(sol.u[i].u,sol.u[i].v,markerstrokewidth=0,markersize=3,color=vortexcolors)
    end
    push!(anims,anim)
end

gif(anims[1])

gif(anims[2])

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
