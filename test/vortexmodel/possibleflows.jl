
# grid
Δx = 0.01
xlim = (-1,1)
ylim = (-1,1)
g = PhysicalGrid(xlim,ylim,Δx);

# circle
R = 0.5
circle = Circle(R,Δx)

# flat plate
c = 1.0
plate = Plate(c,Δx)
T = RigidTransform((0.0,0.0),-π/6)
T(plate);

# vortices
v1 = Vortex(0.75,0.0,1.0);
v2 = Vortex(-0.75,0.0,1.0);
v3 = Vortex(plate.x[1]+3e-2*plate.len*cos(plate.α+π/2),plate.y[1]+3e-2*plate.len*sin(plate.α+π/2),1.0);
v4 = Vortex(plate.x[end]+3e-2*plate.len*cos(plate.α+π/2),plate.y[end]+3e-2*plate.len*sin(plate.α+π/2),1.0);

@testset "Steady unregularized flow around a body" begin
    model = VortexModel(g,bodies=circle)
    ψ = computeψ(model,Ub=(0.0,0.0),U∞=(1.0,0.0),Γb=0.0);
end

@testset "Steady regularized flow around a body" begin
    model = VortexModel(g,bodies=plate,edges=[length(plate)])
    ψ = computeψ(model,U∞=(1.0,0.0));
end

@testset "Unsteady flow without a body" begin
    model = VortexModel(g,vortices=[v1,v2])
    ψ = computeψ(model,U∞=(1.0,1.0));
end

@testset "Unsteady unregularized flow around a body" begin
    model = VortexModel(g,vortices=[v1,v2],bodies=circle)
    ψ = computeψ(model,U∞=(1.0,1.0));
end

@testset "Unsteady regularized flow around a body" begin
    model = VortexModel(g,vortices=[v3,v4],bodies=plate,edges=[1,length(plate)]);
    ψ = computeψ(model,U∞=(1.0,0.0),σ=[SuctionParameterRange(0.0,0.0),SuctionParameterRange(0.0,0.0)]);
end
