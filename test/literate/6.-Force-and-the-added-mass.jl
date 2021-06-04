#=
# 6. Force and the added mass

## Force

This part introduces the calculation of the force and moment on the body through the negative rate of change of impulse in the fluid. The continuous expressions for linear and angular impulse (about the origin) are, in two dimensions,

$\begin{align}
    \boldsymbol{P} &= \int_{V_f} \boldsymbol{x}\times \boldsymbol{\omega}\,\mathrm{d}V + \int_{S_b} \boldsymbol{x}\times \left( \boldsymbol{n}\times\boldsymbol{v}\right)\,\mathrm{d}S \\
    \boldsymbol{\Pi}_0 &= \frac{1}{2}\int_{V_f} \boldsymbol{x}\times \left(\boldsymbol{x}\times\boldsymbol{\omega}\right)\,\mathrm{d}V + \frac{1}{2} \int_{S_b} \boldsymbol{x}\times \left[ \boldsymbol{x}\times\left( \boldsymbol{n} \times \boldsymbol{v} \right) \right] \,\mathrm{d}S
\end{align}$

If there is only a single body, then the force and moment (about the origin) exerted by the fluid on that body are given by

$\boldsymbol{F} = -\rho \frac{\mathrm{d}\boldsymbol{P}}{\mathrm{d}t}, \qquad \boldsymbol{M}_{0} = -\rho \frac{\mathrm{d}\boldsymbol{\Pi}_{0}}{\mathrm{d}t},$

where $\rho$ is the fluid density. In the two-dimensional applications of this package, angular impulse and the moment have only a single component, e.g., $\boldsymbol{\Pi}_0 =\Pi_0\boldsymbol{e}_z$, where $\boldsymbol{e}_z$ is the unit vector out of the plane.
=#

#=
#### Point vortices past a cylinder

To verify our method of calculating the impulse, we create a model of two point vortices of equal and opposite circulation positioned at both sides of the $x$-axis, with the $x$-axis as axis of symmetry, and left of a circular cylinder positioned at the origin. If the top vortex has a positive circulation, the vortices propel each other towards the cylinder and convect past it. This configuration has an analytical solution for the trajecory and impulse.
=#

#md # See the notebook in the examples folder for the analytical solution.

# We use a rectangular grid.
#md # ```@setup 6.-Force-and-the-added-mass
#md # using GridPotentialFlow
#md # using Plots
#md # ```
#!md using GridPotentialFlow
#!md using Plots
Lx = 10.0
Ly = 4.0
xlim = (-Lx/2,Lx/2)
ylim = (-Ly/2,Ly/2)
g = PhysicalGrid(xlim, ylim, 0.04)

#
Rc = 1
circle = Circle(Rc,2*cellsize(g))
pfb = PotentialFlowBody(circle)

Δs = dlength(circle);

# The initial spacing between the vortices is `y∞`. After the vortices pass the cylinder, they should return to this spacing.
y∞ = Rc/2
v1 = Vortex(-Lx/2+10*cellsize(g),y∞,1.0);
v2 = Vortex(-Lx/2+10*cellsize(g),-y∞,-1.0);

# We create the vortex model with these two point vortices and the circle and advance the position of the point vortices over some time. During the time stepping, we compute the impulse with `impulse` and store its history.
model = VortexModel(g,bodies=[pfb],vortices=[v1,v2])
sol = solve(model);

#md # ```@setup 6.-Force-and-the-added-mass
#md # # Analytical trajectory
#md # x(r) = sqrt(r^2-(y∞^2*(r^2-1)^2)/((r^2-1)^2-4*y∞^2))
#md # y(r) = y∞*(r^2-1)/sqrt((r^2-1)^2-4*y∞^2)
#md # s = 0:0.001:π/2
#md # rmin = 1.4375649
#md # c1 = rmin/(-Lx/2+rmin)
#md # c2 = -Lx/2*c1
#md # r = c2 ./ (cos.(s).-c1)
#md # x_trajectory = -reverse(x.(r)); append!(x_trajectory,x.(r))
#md # y_trajectory_upper = reverse(y.(r)); append!(y_trajectory_upper,y.(r))
#md # y_trajectory_lower = -y_trajectory_upper;
#md # # Analytical impulse
#md # Px_func(x,y) = y - y/(x^2+y^2) + y - y/(x^2+y^2);
#md # Py_func(x,y) = -x + x/(x^2+y^2) + x - x/(x^2+y^2);
#md # Px_exact = Px_func.(x_trajectory,y_trajectory_upper);
#md # Py_exact = Py_func.(x_trajectory,y_trajectory_upper);
#md # ```

#!md # Analytical trajectory
#!md x(r) = sqrt(r^2-(y∞^2*(r^2-1)^2)/((r^2-1)^2-4*y∞^2))
#!md y(r) = y∞*(r^2-1)/sqrt((r^2-1)^2-4*y∞^2)
#!md s = 0:0.001:π/2
#!md rmin = 1.4375649
#!md c1 = rmin/(-Lx/2+rmin)
#!md c2 = -Lx/2*c1
#!md r = c2 ./ (cos.(s).-c1)
#!md x_trajectory = -reverse(x.(r)); append!(x_trajectory,x.(r))
#!md y_trajectory_upper = reverse(y.(r)); append!(y_trajectory_upper,y.(r))
#!md y_trajectory_lower = -y_trajectory_upper;
#!md # Analytical impulse
#!md Px_func(x,y) = y - y/(x^2+y^2) + y - y/(x^2+y^2);
#!md Py_func(x,y) = -x + x/(x^2+y^2) + x - x/(x^2+y^2);
#!md Px_exact = Px_func.(x_trajectory,y_trajectory_upper);
#!md Py_exact = Py_func.(x_trajectory,y_trajectory_upper);

Δt = 0.1
T = 0:Δt:64.0
X_hist = []
Px_numerical_hist = []
Py_numerical_hist = []

for t in T
    Ẋ = vortexvelocities!(model)
    X = getvortexpositions(model)
    X = X + Ẋ*Δt
    setvortexpositions!(model,X)

    push!(X_hist,X)

    Px, Py = impulse(model)
    push!(Px_numerical_hist,Px)
    push!(Py_numerical_hist,Py)
end

# When we compare the trajectories and impulse history, the numerical and anaylytical solution should match closely, which is indeed the case.
plot(circle,fillcolor=:black,fillrange=0,fillalpha=0.25,linecolor=:black,linewidth=2,xlabel="x",ylabel="y")
scatter!(model.vortices.x,model.vortices.y,color=:red)
plot!(x_trajectory,y_trajectory_upper,linecolor=:red,label="exact")
plot!(x_trajectory,y_trajectory_lower,linecolor=:red,label="")
plot!((X->X[1]).(X_hist),(X->X[3]).(X_hist),color=:blue,linestyle=:dash,label="simulated")
plot!((X->X[2]).(X_hist),(X->X[4]).(X_hist),color=:blue,linestyle=:dash,label="")

#
plot(x_trajectory,Px_exact,color=:red,label="exact",xlabel="x",ylabel="Px")
plot!((X->X[1]).(X_hist),Px_numerical_hist,color=:blue,linestyle=:dash,label="simulated")

#=
#### Flat plate

We can apply this method also to our flat plate example and compare it to the Biot-Savart method from the `PotentialFlow.jl` package.
=#

#md # See the notebook in the examples folder for the Biot-Savart solution.

# Because we create the `PotentialFlow.jl` model with a moving plate instead of a uniform flow, we specify a translational velocity. The `GridPotentialFlow.jl` model will still use a moving body coordinate system and thus a uniform flow that equals the negative translational velocity.

c = 1.0;   # chord length

ċ = -1.0;  # translational velocity
α = -π/3;   # angle of attack

σLE = 0.0;
σTE = 0.0;

Δt = 5e-2;
tf = 1.0;
T = 0.0:Δt:tf;

#md # ```@setup 6.-Force-and-the-added-mass
#md # # Potential Flow
#md # using PotentialFlow
#md # function compute_ẋ!(ẋ, x, t)
#md #     plate, ambient_sys = x
#md #     motion = ẋ[1]
#md #     motion.ċ, motion.c̈, motion.α̇ = motion.kin(t)
#md #     Plates.enforce_no_flow_through!(plate, motion, ambient_sys, t)
#md #     reset_velocity!(ẋ, x)
#md #     self_induce_velocity!(ẋ, x, t)
#md # end
#md # function shed_new_vorticity!(blobs, plate, motion, t, lesp = 0.0, tesp = 0.0)
#md #     z₊ = (blobs[end-1].z + 2plate.zs[end])/3
#md #     z₋ = (blobs[end].z + 2plate.zs[1])/3
#md #     blob₊ = PotentialFlow.Vortex.Blob(z₊, 1.0, δ)
#md #     blob₋ = PotentialFlow.Vortex.Blob(z₋, 1.0, δ)
#md #     Plates.enforce_no_flow_through!(plate, motion, blobs, t)
#md #     Γ₊, Γ₋, _, _ = Plates.vorticity_flux!(plate, blob₊, blob₋, t, lesp, tesp);
#md #     push!(blobs, PotentialFlow.Vortex.Blob(z₊, Γ₊, blobs[1].δ), PotentialFlow.Vortex.Blob(z₋, Γ₋, blobs[1].δ))
#md # end
#md # N = 128 # number of plate control points (distributed along a extrema Chebyshev grid)
#md # δ = 0.01
#md # plate = PotentialFlow.Plate(N, c, zero(ComplexF64), α)
#md # motion = Plates.RigidBodyMotion(ċ, 0.0);
#md # # Initial step
#md # Δz₀ = im*3Δt*exp(im*plate.α) # vectors perpendicular to the plate
#md # z₋, z₊ = plate.zs[[1,N]] # LE and TE
#md # blobs = PotentialFlow.Vortex.Blob.(Δz₀ .+ [z₊, z₋], 1.0, δ) # First two point vortices are placed close to the LE and TE with unit strength
#md # Plates.enforce_no_flow_through!(plate, motion, (), 0)
#md # Γ₊, Γ₋, _, _ = Plates.vorticity_flux!(plate, blobs[1], blobs[2], 0.0, σLE, σTE); # Determine strength of first two vortices
#md # blobs = PotentialFlow.Vortex.Blob.(Δz₀ .+ [z₊, z₋], [Γ₊, Γ₋], δ) # Create first two point vortices now with calculated strengths
#md # sys₀ = (plate, blobs)
#md # sys = deepcopy(sys₀)
#md # sys₊ = deepcopy(sys₀) # Used for storage during time-marching
#md # ẋs = (motion, allocate_velocity(blobs))
#md # imp = ComplexF64[];
#md # t_hist = Float64[]
#md # global t = 0
#md # push!(t_hist,t)
#md # # Time stepping
#md # for tloc in T[2:end]
#md #     global t += Δt
#md #     global sys
#md #     global sys₊
#md #     push!(imp,Elements.impulse(sys))
#md #     push!(t_hist,t)
#md #     local plate, ambient_ω = sys
#md #     local motion, ambient_u = ẋs
#md #     resize!(sys₊[2], length(sys[2]))
#md #     forward_euler!(sys₊, sys, t, Δt, compute_ẋ!, advect!, ẋs)
#md #     sys, sys₊ = sys₊, sys
#md #     shed_new_vorticity!(sys[2], sys[1], ẋs[1], t, σLE, σTE)
#md # end
#md # force = -diff(imp)/Δt;
#md # ```

#!md # Potential Flow
#!md using PotentialFlow
#!md function compute_ẋ!(ẋ, x, t)
#!md     plate, ambient_sys = x
#!md     motion = ẋ[1]
#!md     motion.ċ, motion.c̈, motion.α̇ = motion.kin(t)
#!md     Plates.enforce_no_flow_through!(plate, motion, ambient_sys, t)
#!md     reset_velocity!(ẋ, x)
#!md     self_induce_velocity!(ẋ, x, t)
#!md end
#!md function shed_new_vorticity!(blobs, plate, motion, t, lesp = 0.0, tesp = 0.0)
#!md     z₊ = (blobs[end-1].z + 2plate.zs[end])/3
#!md     z₋ = (blobs[end].z + 2plate.zs[1])/3
#!md     blob₊ = PotentialFlow.Vortex.Blob(z₊, 1.0, δ)
#!md     blob₋ = PotentialFlow.Vortex.Blob(z₋, 1.0, δ)
#!md     Plates.enforce_no_flow_through!(plate, motion, blobs, t)
#!md     Γ₊, Γ₋, _, _ = Plates.vorticity_flux!(plate, blob₊, blob₋, t, lesp, tesp);
#!md     push!(blobs, PotentialFlow.Vortex.Blob(z₊, Γ₊, blobs[1].δ), PotentialFlow.Vortex.Blob(z₋, Γ₋, blobs[1].δ))
#!md end
#!md N = 128 # number of plate control points (distributed along a extrema Chebyshev grid)
#!md δ = 0.01
#!md plate = PotentialFlow.Plate(N, c, zero(ComplexF64), α)
#!md motion = Plates.RigidBodyMotion(ċ, 0.0);
#!md # Initial step
#!md Δz₀ = im*3Δt*exp(im*plate.α) # vectors perpendicular to the plate
#!md z₋, z₊ = plate.zs[[1,N]] # LE and TE
#!md blobs = PotentialFlow.Vortex.Blob.(Δz₀ .+ [z₊, z₋], 1.0, δ) # First two point vortices are placed close to the LE and TE with unit strength
#!md Plates.enforce_no_flow_through!(plate, motion, (), 0)
#!md Γ₊, Γ₋, _, _ = Plates.vorticity_flux!(plate, blobs[1], blobs[2], 0.0, σLE, σTE); # Determine strength of first two vortices
#!md blobs = PotentialFlow.Vortex.Blob.(Δz₀ .+ [z₊, z₋], [Γ₊, Γ₋], δ) # Create first two point vortices now with calculated strengths
#!md sys₀ = (plate, blobs)
#!md sys = deepcopy(sys₀)
#!md sys₊ = deepcopy(sys₀) # Used for storage during time-marching
#!md ẋs = (motion, allocate_velocity(blobs))
#!md imp = ComplexF64[];
#!md t_hist = Float64[]
#!md global t = 0
#!md push!(t_hist,t)
#!md # Time stepping
#!md for tloc in T[2:end]
#!md     global t += Δt
#!md     global sys
#!md     global sys₊
#!md     push!(imp,Elements.impulse(sys))
#!md     push!(t_hist,t)
#!md     local plate, ambient_ω = sys
#!md     local motion, ambient_u = ẋs
#!md     resize!(sys₊[2], length(sys[2]))
#!md     forward_euler!(sys₊, sys, t, Δt, compute_ẋ!, advect!, ẋs)
#!md     sys, sys₊ = sys₊, sys
#!md     shed_new_vorticity!(sys[2], sys[1], ẋs[1], t, σLE, σTE)
#!md end
#!md force = -diff(imp)/Δt;

Δx = 0.01
xlim = (-0.5,2)
ylim = (-1,1)
g = PhysicalGrid(xlim,ylim,Δx);

#

dsdx = 2
plate = RigidBodyTools.Plate(c,dsdx*Δx)
pfb = PotentialFlowBody(plate,edges=[1,length(plate)],σ=[SuctionParameter(σLE),SuctionParameter(σTE)])
transform = RigidTransform((0.0,0.0),-π/3)
transform(plate);
Δs = dlength(plate)
maximum(Δs/Δx);

# As before, we create the initial vortices based on the time step and uniform flow and define a function that create the new point vortices for the subsequent time steps.

firstvLE = GridPotentialFlow.Vortex(plate.x[1]+3Δt*(-ċ)*cos(plate.α+π/2),plate.y[1]+3Δt*(-ċ)*sin(plate.α+π/2),0.0);
firstvTE = GridPotentialFlow.Vortex(plate.x[end]+3Δt*(-ċ)*cos(plate.α+π/2),plate.y[end]+3Δt*(-ċ)*sin(plate.α+π/2),0.0);

#

function createsheddedvortices(plate,oldvortices)

    vLE = GridPotentialFlow.Vortex(2/3*plate.x[1]+1/3*oldvortices[end-1].x,2/3*plate.y[1]+1/3*oldvortices[end-1].y,0.0)
    vTE = GridPotentialFlow.Vortex(2/3*plate.x[end]+1/3*oldvortices[end].x,2/3*plate.y[end]+1/3*oldvortices[end].y,0.0)

    return vLE,vTE
end

# The model free stream should be the negative of the translational velocity of the plate.

model = VortexModel(g,bodies=[pfb],vortices=[firstvLE,firstvTE],U∞=(-ċ,0.0));

# Note that we have to set the strengths of the new vortices ourselves before calling the `impulse` function. In this case, this is done for us in the `vortexvelocities!`.

Px_hist = Float64[];
Py_hist = Float64[];

# Then we enter a time loop and record the impulse every time we create new vortices.

for tloc in T[2:end]
    Ẋ = vortexvelocities!(model)
    Px, Py = impulse(model)
    X = getvortexpositions(model)
    X = X + Ẋ*Δt
    setvortexpositions!(model,X)
    vLE, vTE = createsheddedvortices(plate,model.vortices)
    pushvortices!(model,vLE,vTE)
    push!(Px_hist,Px)
    push!(Py_hist,Py)
end

# The force history can be obtained by finite differencing the impulse history.

Fx_hist = -diff(Px_hist)/Δt;
Fy_hist = -diff(Py_hist)/Δt;

# We can now compare the positions of the point vortices by shifting the `PotentialFlow.jl` solution by `-tf*ċ` such that the origin of the frame of reference coincides with the center plate. Superimposingt the `GridPotentialFlow.jl` solution then shows that the positions of the point vortices agree very well.

plot(plate,fillcolor=:black,fillrange=0,fillalpha=0.25,linecolor=:black,linewidth=2,xlim=xlim,ylim=ylim,xlabel="x",ylabel="y")
scatter!(real.((v->v.z).(sys[2])).-tf*ċ,imag.((v->v.z).(sys[2])),color=:red,markersize=4,label="PotentialFlow.jl")
scatter!(model.vortices.x,model.vortices.y,color=:blue,markersize=2,label="GridPotentialFlow.jl")

# The vertical impulse and the vertical force (lift) can also be compared and show good agreement as well.
plot(xlabel="t",ylabel="Py")
plot!(T[1:end-1],imag.(imp),color=:blue,label="PotentialFlow.jl")
plot!(T[1:end-1],Py_hist,color=:red,label="GridPotentialFlow.jl")

#
plot(xlabel="t",ylabel="Fy")
plot!(T[2:end-1],imag.(force),color=:blue,label="PotentialFlow.jl")
plot!(T[2:end-1],Fy_hist,color=:red,label="GridPotentialFlow.jl")

#=
## Added mass
=#

# The added mass tensor provides a measure of the inertial influence of the fluid on the body in response to changes in the body's translational or rotational motion. The coefficients of the added mass tensor of a body are obtained by computing the impulse components associated with a unit-valued component of motion. The motion's influence is both direct, via the surface velocity, and indirect, in the bound vortex sheet that develops on the surface.

#=
The added mass for simple geometries can easiliy be calculated in an analytical way. For example, the entries of the translational added mass tensor for an ellipse with semi-major axis $a$ and semi-minor axis $b$ are $m_{xx} = \rho \pi b^2$ for motion in the $x$ direction, $m_{yy} = \rho \pi a^2$ for motion in the $y$-direction, with the diagonal entries $m_{xy}=m_{yx}=0$. By using `addedmass` from this package we can approximate these results numerically. In this case, we have one body and the method will return a matrix with the entries

$\begin{bmatrix}
m_{xx} & m_{xy} \\
m_{yx} & m_{yy}
\end{bmatrix}.$

=#

#md # ```@setup 6.-Force-and-the-added-mass
#md # Lx = 2.0
#md # xlim = (-Lx/2,Lx/2)
#md # ylim = (-Lx/2,Lx/2)
#md # g = PhysicalGrid(xlim,ylim,0.01)
#md # Δx = cellsize(g);
#md # ```
#!md Lx = 2.0
#!md xlim = (-Lx/2,Lx/2)
#!md ylim = (-Lx/2,Lx/2)
#!md g = PhysicalGrid(xlim,ylim,0.01)
#!md Δx = cellsize(g);


a = 0.5
b = 0.25
ellipse = Ellipse(a,b,Δx)
pfb = PotentialFlowBody(ellipse)
model = VortexModel(g,bodies=[pfb])
M = GridPotentialFlow.addedmass(model)

# We can compare the values using the `@test` macro.
using Test
@test isapprox(M[1,1], π*b^2, rtol=1e-1)
#
@test isapprox(M[2,2], π*a^2, atol=1e-1)

# In the case of $N$ bodies, `addedmass` returns the translational added mass tensor of size $2N$-by-$2N$. As an example, we will compute the added mass tensor for an array of cylinders. To compare it with an analytically obtained solution, we will calculate the added mass matrix tensor for the following gap-to-radius ratios.

GRratios = [0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.2,1.4,1.6];

#md # ```@setup 6.-Force-and-the-added-mass
#md # function rectangulararray(rows,columns,spacing)
#md #     centers = VectorData(rows*columns)
#md #     for j in 1:columns
#md #         for i in 1:rows
#md #             xc = (i-1)*spacing
#md #             yc = (j-1)*spacing
#md #             centers.u[(j-1)*rows+i] = xc
#md #             centers.v[(j-1)*rows+i] = yc
#md #         end
#md #     end
#md #     centers.u .-= mean(centers.u)
#md #     centers.v .-= mean(centers.v)
#md #     return centers
#md # end
#md # ```

#!md function rectangulararray(rows,columns,spacing)
#!md      centers = VectorData(rows*columns)
#!md      for j in 1:columns
#!md         for i in 1:rows
#!md             xc = (i-1)*spacing
#!md             yc = (j-1)*spacing
#!md             centers.u[(j-1)*rows+i] = xc
#!md             centers.v[(j-1)*rows+i] = yc
#!md         end
#!md     end
#!md     centers.u .-= mean(centers.u)
#!md     centers.v .-= mean(centers.v)
#!md     return centers
#!md  end

# For this case, the array consists of three rows of three cylinders.

n = 100
rows = 3
columns = 3
N = rows*columns
R = 0.20
𝒱 = π*R^2
bodies = fill(Circle(R,2*cellsize(g)),N);
Δs = minimum(dlength(bodies[1]));

# We loop over the gap-to-radius ratios, position the bodies, and compute the ratio of the largest eigenvalue to the largest diagonal element of the translational added mass tensor. To position the bodies, we defined the method `rectangulararray` (see notebook in the examples folder).
using Statistics: mean
using LinearAlgebra: eigen
#
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

# The same ratios, obtained in an analytical way, are available in literature [^1] and can be used to verify the numerical values.
nineRodArrayChen = [2.3637,2.2092,2.1007,2.0120,1.9350,1.8665,1.7494,1.6531,1.5732,1.5066,1.4508];
plot(xlabel="GR",ylabel="λ/max(Mij)")
plot!(GRratios,nineRodArrayChen,label="Chen1975")
plot!(GRratios,λoverMratios,label="GridPotentialFlow.jl")

# [^1]: Chen, S. S. (1975) "Vibration of nuclear fuel bundles," *Nuclear Engineering and Design*, 35 (3), 399-–422.
