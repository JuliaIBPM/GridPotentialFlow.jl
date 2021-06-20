using GridPotentialFlow
using Plots
Lx = 10.0
Ly = 4.0
xlim = (-Lx/2,Lx/2)
ylim = (-Ly/2,Ly/2)
g = PhysicalGrid(xlim, ylim, 0.04)

Rc = 1
circle = Circle(Rc,2*cellsize(g))
pfb = PotentialFlowBody(circle)

Δs = dlength(circle);

y∞ = Rc/2
v1 = Vortex(-Lx/2+10*cellsize(g),y∞,1.0);
v2 = Vortex(-Lx/2+10*cellsize(g),-y∞,-1.0);

x(r) = sqrt(r^2-(y∞^2*(r^2-1)^2)/((r^2-1)^2-4*y∞^2))
y(r) = y∞*(r^2-1)/sqrt((r^2-1)^2-4*y∞^2)
s = 0:0.001:π/2
rmin = 1.4375649
c1 = rmin/(-Lx/2+rmin)
c2 = -Lx/2*c1
r = c2 ./ (cos.(s).-c1)
x_trajectory = -reverse(x.(r)); append!(x_trajectory,x.(r))
y_trajectory_upper = reverse(y.(r)); append!(y_trajectory_upper,y.(r))
y_trajectory_lower = -y_trajectory_upper;

Px_func(x,y) = y - y/(x^2+y^2) + y - y/(x^2+y^2);
Py_func(x,y) = -x + x/(x^2+y^2) + x - x/(x^2+y^2);
Px_exact = Px_func.(x_trajectory,y_trajectory_upper);
Py_exact = Py_func.(x_trajectory,y_trajectory_upper);

model = VortexModel(g,bodies=[pfb],vortices=[v1,v2]);

function rhs(X,model,t)
    setvortexpositions!(model,X)
    Ẋ = vortexvelocities!(model)
    return Ẋ
end

import OrdinaryDiffEq
X = getvortexpositions(model)
prob = OrdinaryDiffEq.ODEProblem(rhs,X,(0.0,64.0),model);
sol = OrdinaryDiffEq.solve(prob,dt=0.1,OrdinaryDiffEq.RK4(),dense=false,adaptive=false);

Px_numerical_hist = []
Py_numerical_hist = []
for u in sol.u
    setvortexpositions!(model,u)
    Px, Py = impulse(model)
    push!(Px_numerical_hist,Px)
    push!(Py_numerical_hist,Py)
end

plot(circle,fillcolor=:black,fillrange=0,fillalpha=0.25,linecolor=:black,linewidth=2,xlabel="x",ylabel="y")
scatter!(model.vortices.x,model.vortices.y,color=:red)
plot!(x_trajectory,y_trajectory_upper,linecolor=:red,label="exact")
plot!(x_trajectory,y_trajectory_lower,linecolor=:red,label="")
plot!(map(s->s.u[1],sol.u),map(s->s.v[1],sol.u),color=:blue,linestyle=:dash,label="simulated")
plot!(map(s->s.u[2],sol.u),map(s->s.v[2],sol.u),color=:blue,linestyle=:dash,label="")

plot(x_trajectory,Px_exact,color=:red,label="exact",xlabel="x",ylabel="Px")
plot!(map(s->s.u[1],sol.u),Px_numerical_hist,color=:blue,linestyle=:dash,label="simulated")

c = 1.0;   # chord length

ċ = -1.0;  # translational velocity
α = -π/3;   # angle of attack

σLE = 0.0;
σTE = 0.0;

Δt = 5e-2;
tf = 1.0;
T = 0.0:Δt:tf;

using PotentialFlow
function compute_ẋ!(ẋ, x, t)
    plate, ambient_sys = x
    motion = ẋ[1]
    motion.ċ, motion.c̈, motion.α̇ = motion.kin(t)
    Plates.enforce_no_flow_through!(plate, motion, ambient_sys, t)
    reset_velocity!(ẋ, x)
    self_induce_velocity!(ẋ, x, t)
end
function shed_new_vorticity!(blobs, plate, motion, t, lesp = 0.0, tesp = 0.0)
    z₊ = (blobs[end-1].z + 2plate.zs[end])/3
    z₋ = (blobs[end].z + 2plate.zs[1])/3
    blob₊ = PotentialFlow.Vortex.Blob(z₊, 1.0, δ)
    blob₋ = PotentialFlow.Vortex.Blob(z₋, 1.0, δ)
    Plates.enforce_no_flow_through!(plate, motion, blobs, t)
    Γ₊, Γ₋, _, _ = Plates.vorticity_flux!(plate, blob₊, blob₋, t, lesp, tesp);
    push!(blobs, PotentialFlow.Vortex.Blob(z₊, Γ₊, blobs[1].δ), PotentialFlow.Vortex.Blob(z₋, Γ₋, blobs[1].δ))
end
N = 128 # number of plate control points (distributed along a extrema Chebyshev grid)
δ = 0.01
plate = PotentialFlow.Plate(N, c, zero(ComplexF64), α)
motion = Plates.RigidBodyMotion(ċ, 0.0);

Δz₀ = im*3Δt*exp(im*plate.α) # vectors perpendicular to the plate
z₋, z₊ = plate.zs[[1,N]] # LE and TE
blobs = PotentialFlow.Vortex.Blob.(Δz₀ .+ [z₊, z₋], 1.0, δ) # First two point vortices are placed close to the LE and TE with unit strength
Plates.enforce_no_flow_through!(plate, motion, (), 0)
Γ₊, Γ₋, _, _ = Plates.vorticity_flux!(plate, blobs[1], blobs[2], 0.0, σLE, σTE); # Determine strength of first two vortices
blobs = PotentialFlow.Vortex.Blob.(Δz₀ .+ [z₊, z₋], [Γ₊, Γ₋], δ) # Create first two point vortices now with calculated strengths
sys₀ = (plate, blobs)
sys = deepcopy(sys₀)
sys₊ = deepcopy(sys₀) # Used for storage during time-marching
ẋs = (motion, allocate_velocity(blobs))
imp = ComplexF64[];
t_hist = Float64[]
global t = 0
push!(t_hist,t)

for tloc in T[2:end]
    global t += Δt
    global sys
    global sys₊
    push!(imp,Elements.impulse(sys))
    push!(t_hist,t)
    local plate, ambient_ω = sys
    local motion, ambient_u = ẋs
    resize!(sys₊[2], length(sys[2]))
    forward_euler!(sys₊, sys, t, Δt, compute_ẋ!, advect!, ẋs)
    sys, sys₊ = sys₊, sys
    shed_new_vorticity!(sys[2], sys[1], ẋs[1], t, σLE, σTE)
end
force = -diff(imp)/Δt;

Δx = 0.01
xlim = (-0.5,2)
ylim = (-1,1)
g = PhysicalGrid(xlim,ylim,Δx);

dsdx = 2
plate = RigidBodyTools.Plate(c,dsdx*Δx)
pfb = PotentialFlowBody(plate,edges=[1,length(plate)],σ=[SuctionParameter(σLE),SuctionParameter(σTE)])
transform = RigidTransform((0.0,0.0),-π/3)
transform(plate);
Δs = dlength(plate)
maximum(Δs/Δx);

firstvLE = GridPotentialFlow.Vortex(plate.x[1]+3Δt*(-ċ)*cos(plate.α+π/2),plate.y[1]+3Δt*(-ċ)*sin(plate.α+π/2),0.0);
firstvTE = GridPotentialFlow.Vortex(plate.x[end]+3Δt*(-ċ)*cos(plate.α+π/2),plate.y[end]+3Δt*(-ċ)*sin(plate.α+π/2),0.0);

function createsheddedvortices(plate,oldvortices)

    vLE = GridPotentialFlow.Vortex(2/3*plate.x[1]+1/3*oldvortices[end-1].x,2/3*plate.y[1]+1/3*oldvortices[end-1].y,0.0)
    vTE = GridPotentialFlow.Vortex(2/3*plate.x[end]+1/3*oldvortices[end].x,2/3*plate.y[end]+1/3*oldvortices[end].y,0.0)

    return vLE,vTE
end

model = VortexModel(g,bodies=[pfb],vortices=[firstvLE,firstvTE],U∞=(-ċ,0.0));

Px_hist = Float64[];
Py_hist = Float64[];
sol = solve(model);
for tloc in T[2:end]
    X = getvortexpositions(model) # gets bigger every time step because we add vortices
    Ẋ = deepcopy(X)

    solve!(sol, model)
    setvortexstrengths!(model, sol.δΓ_vec, length(X.u)-1:length(X.u))
    subtractcirculation!(model.bodies, sol.δΓ_vec)
    Px, Py = impulse(model)

    vortexvelocities!(Ẋ, model, sol.ψ)
    X .= X .+ Ẋ*Δt
    setvortexpositions!(model, X)

    vLE, vTE = createsheddedvortices(plate,model.vortices)
    pushvortices!(model,vLE,vTE)
    push!(Px_hist,Px)
    push!(Py_hist,Py)
end

Fx_hist = -diff(Px_hist)/Δt;
Fy_hist = -diff(Py_hist)/Δt;

plot(plate,fillcolor=:black,fillrange=0,fillalpha=0.25,linecolor=:black,linewidth=2,xlim=xlim,ylim=ylim,xlabel="x",ylabel="y")
scatter!(real.((v->v.z).(sys[2])).-tf*ċ,imag.((v->v.z).(sys[2])),color=:red,markersize=4,label="PotentialFlow.jl")
scatter!(model.vortices.x,model.vortices.y,color=:blue,markersize=2,label="GridPotentialFlow.jl")

plot(xlabel="t",ylabel="Py")
plot!(T[1:end-1],imag.(imp),color=:blue,label="PotentialFlow.jl")
plot!(T[1:end-1],Py_hist,color=:red,label="GridPotentialFlow.jl")

plot(xlabel="t",ylabel="Fy")
plot!(T[2:end-1],imag.(force),color=:blue,label="PotentialFlow.jl")
plot!(T[2:end-1],Fy_hist,color=:red,label="GridPotentialFlow.jl")

Lx = 2.0
xlim = (-Lx/2,Lx/2)
ylim = (-Lx/2,Lx/2)
g = PhysicalGrid(xlim,ylim,0.01)
Δx = cellsize(g);


a = 0.5
b = 0.25
ellipse = Ellipse(a,b,Δx)
pfb = PotentialFlowBody(ellipse)
model = VortexModel(g,bodies=[pfb])
M = GridPotentialFlow.addedmass(model)

using Test
@test isapprox(M[1,1], π*b^2, rtol=1e-1)

@test isapprox(M[2,2], π*a^2, atol=1e-1)

GRratios = [0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.2,1.4,1.6];


function rectangulararray(rows,columns,spacing)
     centers = VectorData(rows*columns)
     for j in 1:columns
        for i in 1:rows
            xc = (i-1)*spacing
            yc = (j-1)*spacing
            centers.u[(j-1)*rows+i] = xc
            centers.v[(j-1)*rows+i] = yc
        end
    end
    centers.u .-= mean(centers.u)
    centers.v .-= mean(centers.v)
    return centers
 end

n = 100
rows = 3
columns = 3
N = rows*columns
R = 0.20
𝒱 = π*R^2
bodies = fill(Circle(R,2*cellsize(g)),N);
Δs = minimum(dlength(bodies[1]));

using Statistics: mean
using LinearAlgebra: eigen

λoverMratios = zeros(length(GRratios))
for idx in 1:length(GRratios)
    global bodies
    gap = GRratios[idx]*R
    spacing = 2*R+gap
    bodycenters = rectangulararray(rows,columns,spacing)
    for i in 1:N
        Tf = RigidTransform((bodycenters.u[i],bodycenters.v[i]),0.0)
        global bodies[i] = Tf(deepcopy(bodies[i]))
    end
    body_list = PotentialFlowBody.(bodies)
    vm = VortexModel(g,bodies=body_list);
    M = GridPotentialFlow.addedmass(vm)/𝒱
    eigM = eigen(M);
    max_eig_value_coef = maximum(real(eigM.values))
    max_self_added_mass_coef = maximum(M)
    λoverMratios[idx] = max_eig_value_coef/max_self_added_mass_coef
end

nineRodArrayChen = [2.3637,2.2092,2.1007,2.0120,1.9350,1.8665,1.7494,1.6531,1.5732,1.5066,1.4508];
plot(xlabel="GR",ylabel="λ/max(Mij)")
plot!(GRratios,nineRodArrayChen,label="Chen1975")
plot!(GRratios,λoverMratios,label="GridPotentialFlow.jl")

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

