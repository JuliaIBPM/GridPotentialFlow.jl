#=
# 1. Basic potential flow problem

This part introduces how to solve the basic potential flow problem in this package. We consider here a staggered, Cartesian grid with uniform cell size and of infinite extent. The basic (unbounded) potential flow problem is expressed as

$\mathsf{Ls} = -\mathsf{w},$

where $\mathsf{L}$ is the discrete 5-point Laplacian operator, $\mathsf{s}$ is the discrete streamfunction, and $\mathsf{w}$ is the discrete vorticity. Following the vortex-in-cell approach, the discrete vorticity is obtained by regularizing the vorticity from the $N_v$ vortex elements onto the cell vertices,

$\mathsf{w}_{i, j}=\sum_{q=1}^{N_{v}} \frac{1}{\Delta x} \Gamma_{v, q} d\left(\frac{\mathsf{x}_{i}-X_{q}}{\Delta x}\right) d\left(\frac{\mathsf{y}_{j}-Y_{q}}{\Delta x}\right)$

where $d$ is the $M_{4}'$ interpolation kernel, and $(X_{q},Y_{q})$ and $\Gamma_{v,q}$ are the position and strength of the $q$th vortex element.
=#


#=
## Flow around a point vortex
Now we will solve this discrete potential flow problem using `GridPotentialFlow.jl` to obtain the streamfunction field and velocity field around a single point vortex.
=#

# The first step in creating *any* `GridPotentialFlow` model is to create a `PhysicalGrid` to discretize the domain.
using GridPotentialFlow
Δx = 0.01
Lx = 2.0
xlim = (-Lx/2,Lx/2)
ylim = (-Lx/2,Lx/2)
g = PhysicalGrid(xlim,ylim,Δx,optimize=false);

# The second step is to create our point vortex using `Vortex`.
v = Vortex(0.0,0.0,1.0);

# Now we can create a `VortexModel` using the grid and an array containing the point vortex.
model = VortexModel(g,vortices=[v]);

# The discrete streamfunction `s` is then obtained using `streamfunction`.
s = streamfunction(model);
using Plots
plot(s,g,xlabel="x",ylabel="y")

# The function `streamfunction` returns a `Nodes` array. If we want to perform differential calculus operations this data, we can use the methods of the  `CartesianGrids` package. For example, we can easily obtain the velocity field from the streamfunction field using the `curl` operation.
q = curl(s);
plot(q,g,xlabel="x",ylabel="y")

#=
## Accuracy of the discretized Poisson equation
To verify that the discretization technique of the Poisson equation is second-order accurate, we perform a mesh refinement analysis for a flow consisting of point vortices of random strenght that are randomly positioned on the lower-left quadrant of a square domain.
=#

# The error is calculated by comparing the numerically calculated streamfunction field with the exact solution for the streamfunction field. We therefore need a function to calculate the exact streamfunction.
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

# We create four vortices with random strenghts, randomly positioned in the lower-left quadrant of the domain.
nv = 4;
vl = Vortex.(-Lx/4 .+ 0.4*Lx*(rand(nv).-0.5),-Lx/4 .+ 0.4*Lx*(rand(nv).-0.5),0.5*rand(nv).+0.5);

# Next, we create a series of grids, with each grid doubling the number of grid points of the previous grid in each direction.
grids = [PhysicalGrid(xlim,ylim,Lx/(nx-2),optimize=false) for nx in [2^p for p in 5:9]];

# The error is calculated as $\epsilon = \Vert \psi(\mathsf{x},\mathsf{y})/\Delta x - \mathsf{s} \Vert_2 / \Vert \mathsf{s} \Vert_2$, for which we use the `norm` function.
using LinearAlgebra: norm

# We now loop over the grids and calculate the error.
errors = []
gridspacings = []
for g in grids
    ## LGF
    local model = VortexModel(g,vortices=vl)
    s_lgf = streamfunction(model);

    ## Exact solution
    s_exact = Nodes(Dual,size(g))
    s_temp = Nodes(Dual,size(g))
    for v in vl
        ψ_vortex!(s_temp,v,g)
        s_exact += s_temp
    end

    ## Bring to same reference level. A constant value can be added/subtracted to any potential flow solution
    s_lgf .-= s_lgf[g.I0[1],2];
    s_exact .-= s_exact[g.I0[1],2];

    error = s_lgf-s_exact
    idx = g.I0[1] + 1 : g.N[1] - 1 # Only look at top right corner
    push!(errors,norm(error[idx,idx])/norm(s_lgf[idx,idx]))
    push!(gridspacings,g.Δx)
end

# And finally, we create a log-log plot of the error versus the grid spacing.
firstorder=(errors[end]/gridspacings[end].*gridspacings).^1
secondorder=(sqrt(errors[end])/gridspacings[end].*gridspacings).^2
p=plot(gridspacings,errors,xaxis=:log,yaxis=:log,marker=:circle,lab="error",xlabel="dx",legend=:bottomright,title="Error")
plot!(gridspacings,firstorder,lab="1st order",linestyle=:dot,linecolor=:black)
plot!(gridspacings,secondorder,lab="2nd order",linestyle=:dash,linecolor=:black)

#=
## Corotating point vortices
Now we will try advancing a vortex model in time. The simplest unsteady vortex model consists of two point vortices. When both point vortices have the same strength, they will rotate around each other on a trajectory that is easy to describe analytically. In this example, we will compare the analytical and simulated trajectories during one revolution.
=#

# We create two vortices at a distance $d$ from each other and give them a strength $\Gamma$
d = Lx/2
Γ = 1
v1 = Vortex(d/2,0.0,Γ);
v2 = Vortex(-d/2,0.0,Γ);

# We can analytically determine the time $T$ it takes for these vortices to complete one revolution around their centroid.
Vθ = Γ/(2*π*d); # Analytical tangential velocity
T = π*d/Vθ # Analytical period

# Let's now create a vortex model with the two point vortices.
model = VortexModel(g,vortices=[v1,v2]);

# If we update the positions of the point vortices repeatedly to simulate the model advancing in time from $t=0$ to $t=T$, we can check if they end up again at their original positions.
tspan = (0.0,T);

# To step in time, we can simply update the position of the $q$th vortex as $X^{n+1}_q = X^{n}_q + Δt Ẋ^{n}_q$ (forward Euler) in a for-loop. Alternatively, we can make use of the `OrdinaryDiffEq.jl` package and use, for example, their fourth-order Runge-Kutta time stepping scheme.

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

# We finally create the problem and solve it with a fourth-order Runge-Kutta time stepping scheme.
X = getvortexpositions(model)
prob = OrdinaryDiffEq.ODEProblem(rhs,X,tspan,model);
sol = OrdinaryDiffEq.solve(prob,dt=0.1,OrdinaryDiffEq.RK4(),dense=false,adaptive=false);

# When we plot the x- and y-position versus time of the first vortex, we see it completed one revolution.
plot(sol.t,map(s->s.u[1],sol.u),xlabel='t',label="x₁")
plot!(sol.t,map(s->s.v[1],sol.u),label="y₁")

#jl @testset "Corotating point vortices" begin
#jl     @test isapprox(sol.u[end].u[1], 0.5; atol = 1e-2)
#jl     @test isapprox(sol.u[end].u[2], -0.5; atol = 1e-2)
#jl     @test isapprox(sol.u[end].v[1], 0.0; atol = 1e-2)
#jl     @test isapprox(sol.u[end].v[2], 0.0; atol = 1e-2)
#jl end

#=
## Corotating vortex patches
A more complex example is the evolution of two circular regions of spatially uniform vorticity that have an equal radius $r_0$ and circulation $\Gamma$, and whose centers are separated by a distance $d_0$. If $r_0/d_0$ is small enough, these two vortex patches can be regarded as two point vortices with the same circulation and separation. We therefore assume the period of this system is approximately the same as the the period for the point vortices case.
=#

#md # ```@setup 1.-Basic-potential-flow-problem
#md # xlim = (-2,2);
#md # ylim = (-2,2);
#md # Δx = 0.05;
#md # g = PhysicalGrid(xlim,ylim,Δx,optimize=false);
#md # ```
#!md xlim = (-2,2);
#!md ylim = (-2,2);
#!md Δx = 0.05;
#!md g = PhysicalGrid(xlim,ylim,Δx,optimize=false);

# We discretize the vortex patches with point vortices arranged on concentric rings using the following function.

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

# We will simulate two cases with different values for $r_0/d_0$, the ratio of the vortex patch radius to the distance between them.
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

# If the ratio $r_0/d_0$ is big enough, the two vortex patches merge with eachother, a result that has been widely reported in literature [^1].

gif(anims[2])

# [^1]: Eldredge J. D. (2019) "Mathematical  Modeling  of  Unsteady  Inviscid  Flows, Interdisciplinary Applied Mathematics" *Springer*, vol. 50.
