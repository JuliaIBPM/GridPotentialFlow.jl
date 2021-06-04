using GridPotentialFlow
using Plots
Δx = 0.01
Lx = 2.0
xlim = (-Lx/2,Lx/2)
ylim = (-Lx/2,Lx/2)
g = PhysicalGrid(xlim,ylim,Δx);

R = Lx/8;
b_left = PotentialFlowBody(Circle(R,Δx),Γ=1.0)
T = RigidTransform((-Lx/4,0.0),0.0)
T(b_left);
b_right = PotentialFlowBody(Circle(R,Δx))
T = RigidTransform((Lx/4,0.0),0.0)
T(b_right);

model = VortexModel(g,bodies=[b_left,b_right])
s = streamfunction(model)
plot(s,g)
plot!(b_left,fillcolor=:black,fillrange=0,fillalpha=0.25,linecolor=:black,linewidth=2)
plot!(b_right,fillcolor=:black,fillrange=0,fillalpha=0.25,linecolor=:black,linewidth=2)

Δx = 0.004
Lx = 2.0
xlim = (-Lx/2,Lx/2)
ylim = (-Lx/2,Lx/2)
g = PhysicalGrid(xlim,ylim,Δx);

a = Lx/4 # semi-major axis
AR_list = [1.0,2.0,3.0]
ellipses = [Ellipse(a,a/AR,Δx) for AR in AR_list]
bodies = PotentialFlowBody[]
push!(bodies,[PotentialFlowBody(SplinedBody(hcat(e.x,e.y),3*cellsize(g))) for e in ellipses]...)
push!(bodies,PotentialFlowBody(Plate(2*a,3*Δx)));
colors = [:blue,:red,:green,:black]
plot()
for i in 1:length(bodies)
    plot!(bodies[i],fillrange=0,fillalpha=0.0,linecolor=colors[i],linewidth=2)
end
plot!(xlim=xlim,ylim=ylim,xlabel="x",ylabel="y")

f₀_list = []
for i in 1:length(bodies)
    Δs = dlength(bodies[i]);
    model = VortexModel(g,bodies=[bodies[i]]);
    ones = ScalarData(length(bodies[i]))
    ones .= 1.0
    f₀ = model.system.ibp.Sfact\ones
    push!(f₀_list,f₀)
end

append!(f₀_list[end],reverse!(f₀_list[end]));
Δs_plate = dlengthmid(bodies[end]);
append!(Δs_plate,reverse!(Δs_plate));
plot(ylim=(-2,0),xlabel="p/N",title="f₀ for ellipses of different aspect ratios")
for i in 1:length(bodies)-1
    plot!((1:length(f₀_list[i]))/length(f₀_list[i]),f₀_list[i]./dlength(bodies[i]),linecolor=colors[i],label="AR=$(AR_list[i])")
end
plot!((1:length(f₀_list[end]))/length(f₀_list[end]),0.5*f₀_list[end]./Δs_plate,linecolor=colors[end],label="AR=inf")

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

