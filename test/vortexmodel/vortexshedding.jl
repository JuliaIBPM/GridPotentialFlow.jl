
# grid
Δx = 0.005
xlim = (-0.5,1.0)
ylim = (-0.7,0.7)
g = PhysicalGrid(xlim,ylim,Δx)

# flat plate
c = 1.0
plate = Plate(c,100)
T = RigidTransform((0.0,0.0),-π/3)
T(plate);

function createsheddedvortices(plate,oldvortices,Δt)
    vLE = Vortex(2/3*plate.x[1]+1/3*oldvortices[end-1].x,2/3*plate.y[1]+1/3*oldvortices[end-1].y,0.0)
    vTE = Vortex(2/3*plate.x[end]+1/3*oldvortices[end].x,2/3*plate.y[end]+1/3*oldvortices[end].y,0.0)
    return vLE, vTE
end

@testset "Basic vortex shedding" begin

    Δt = 2e-2
    T = 0:Δt:0.2

    vLE = Vortex(plate.x[1]+3Δt*plate.len*cos(plate.α+π/2),plate.y[1]+3Δt*plate.len*sin(plate.α+π/2),1.0)
    vTE = Vortex(plate.x[end]+3Δt*plate.len*cos(plate.α+π/2),plate.y[end]+3Δt*plate.len*sin(plate.α+π/2),1.0)

    model = VortexModel(g,bodies=plate,edges=[1,length(plate)],vortices=[vLE,vTE]);

    for t in T
        Ẋ = computevortexvelocities(model,parameters=ModelParameters(U∞=(1.0,0.0)))
        vortices = deepcopy(model.vortices.list)
        updateposition!.(vortices,Ẋ.u,Ẋ.v,Δt)
        setvortices!(model,vortices)
        vLE, vTE = createsheddedvortices(plate,model.vortices.list,Δt)
        pushvortices!(model,vLE,vTE)
    end
end
