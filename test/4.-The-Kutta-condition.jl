using GridPotentialFlow
using Plots
Δx = 0.01
Lx = 2.0
xlim = (-Lx/2,Lx/2)
ylim = (-Lx/2,Lx/2)
g = PhysicalGrid(xlim,ylim,Δx);

c = Lx/2 # chord length
α = -π/6 # angle of attack
plate = Plate(c,4*cellsize(g))
Tr = RigidTransform((0.0,0.0),α)
Tr(plate)
pfb = PotentialFlowBody(plate);
model = VortexModel(g,bodies=[pfb]);

setU∞(model,(1.0,0.0))
sol = solve(model);

plot(sol.ψ,g,levels=60);
plot!(plate,linecolor=:black,linewidth=2,xlabel="x",ylabel="y")

plot(sol.f,xlabel="body point index",ylabel="f",legend=false)

ones = ScalarData(length(plate))
ones .= 1.0
f₀ = model.system.ibp.Sfact\ones
plot(f₀)

plot(sol.f./f₀,xlabel="body point index",ylabel="f̃",legend=false)

pfb = PotentialFlowBody(plate,edges=[length(plate)])
model = VortexModel(g,bodies=[pfb],U∞=(1.0,0.0))
sol = solve(model);

plot(sol.ψ,g);
plot!(plate,linecolor=:black,linewidth=2,xlabel="x",ylabel="y")

plot(plot(sol.f,xlabel="body point index",ylabel="f",legend=false),plot(sol.f./f₀,xlabel="body point index",ylabel="f̃",legend=false),size=[800,300])

Δt = 1e-2
vTE = Vortex(plate.x[end]+3*Δt*cos(α+π/2),plate.y[end]+3*Δt*sin(α+π/2),0.0);

pfb = PotentialFlowBody(plate,edges=[length(plate)])
model = VortexModel(g,bodies=[pfb],vortices=[vTE],U∞=(1.0,0.0))
sol = solve(model);
plot(sol.ψ,g);
plot!(plate,linecolor=:black,linewidth=2)
scatter!(model.vortices.x,model.vortices.y,color=:black,markersize=2,xlabel="x",ylabel="y")

plot(plot(sol.f,xlabel="body point index",ylabel="f",legend=false),plot(sol.f./f₀,xlabel="body point index",ylabel="f̃",legend=false),size=[800,300])

vLE = Vortex(plate.x[1]+3*Δt*plate.len*cos(plate.α+π/2),plate.y[1]+3*Δt*plate.len*sin(plate.α+π/2),0.0);
vTE = Vortex(plate.x[end]+3*Δt*cos(α+π/2),plate.y[end]+3*Δt*sin(α+π/2),0.0);

pfb = PotentialFlowBody(plate,edges=[1,length(plate)])
model = VortexModel(g,bodies=[pfb],vortices=[vLE,vTE],U∞=(1.0,0.0))
sol = solve(model);
plot(sol.ψ,g);
plot!(plate,linecolor=:black,linewidth=2)
scatter!(model.vortices.x,model.vortices.y,markersize=2,xlabel="x",ylabel="y")

plot(plot(sol.f,xlabel="body point index",ylabel="f",legend=false),plot(sol.f./f₀,xlabel="body point index",ylabel="f̃",legend=false),size=[800,300])

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

