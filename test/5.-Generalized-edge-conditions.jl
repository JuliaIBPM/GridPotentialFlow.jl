using GridPotentialFlow
using Plots
Δx = 0.01
Lx = 2.0
xlim = (-Lx/2,Lx/2)
ylim = (-Lx/2,Lx/2)
g = PhysicalGrid(xlim,ylim,Δx);
c = Lx/2 # chord length
α = -π/3 # angle of attack
plate = Plate(c,4*cellsize(g));
Tr = RigidTransform((0.0,0.0),α)
Tr(plate)
Δs = dlengthmid(plate);

σLE_list = [0.0,0.05,0.1];
pfb_list = [PotentialFlowBody(plate,edges=[1,length(plate)],σ=[SuctionParameter(σLE),SuctionParameter(0.0)]) for σLE in σLE_list];

Δt = 2e-2
vLE = Vortex(plate.x[1]+3Δt*plate.len*cos(plate.α+π/2),plate.y[1]+3Δt*plate.len*sin(plate.α+π/2),0.0);
vTE = Vortex(plate.x[end]+3Δt*cos(α+π/2),plate.y[end]+3Δt*sin(α+π/2),0.0);

function createsheddedvortices(plate,oldvortices)

    vLE = Vortex(2/3*plate.x[1]+1/3*oldvortices[end-1].x,2/3*plate.y[1]+1/3*oldvortices[end-1].y,0.0)
    vTE = Vortex(2/3*plate.x[end]+1/3*oldvortices[end].x,2/3*plate.y[end]+1/3*oldvortices[end].y,0.0)

    return vLE, vTE
end

model_list = [VortexModel(g,bodies=[pfb],vortices=[vLE,vTE],U∞=(1.0,0.0)) for pfb in pfb_list];

T = 0:Δt:0.2
for t in T
    for vm in model_list
        X = getvortexpositions(vm)
        Ẋ = vortexvelocities!(vm)
        X .= X .+ Ẋ*Δt
        setvortexpositions!(vm, X)
        vLEnew, vTEnew = createsheddedvortices(plate,vm.vortices[end-1:end])
        pushvortices!(vm,vLEnew,vTEnew)
    end
end

colors = [:red,:blue,:green];
plot(plate,fillcolor=:black,fillrange=0,fillalpha=0.25,linecolor=:black,linewidth=2,xlabel="x",ylabel="y")
for i in 1:length(model_list)
    plot!(model_list[i].vortices.x[4:2:end],model_list[i].vortices.y[4:2:end],color=colors[i],marker=:circle,markersize=2)
    plot!(model_list[i].vortices.x[3:2:end],model_list[i].vortices.y[3:2:end],color=colors[i],marker=:circle,markersize=2)
    scatter!(model_list[i].vortices.x[1:2],model_list[i].vortices.y[1:2],color=colors[i],marker=:circle,markersize=2,label="σLE=$(σLE_list[i])")
end
plot!()

plot(xlabel="body point index",ylabel="f̃")
for i in 1:length(model_list)
    sol = solve(model_list[i]);
    plot!(sol.f./model_list[i].system.f₀,linecolor=colors[i],label="σLE=$(σLE_list[i])")
end
plot!()

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

