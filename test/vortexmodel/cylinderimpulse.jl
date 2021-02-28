
# grid
Δx = 0.01
xlim = (-1,1)
ylim = (-1,1)
g = PhysicalGrid(xlim,ylim,Δx);

# circle
Rc = 0.5
circle = Circle(Rc,Δx)

# ellipse
a = 0.5
b = 0.25
ellipse = Ellipse(a,b,Δx)

@testset "Impulse" begin
    model = VortexModel(g,bodies=circle)
    Px, Py = computeimpulse(model,parameters=ModelParameters(Ub=(1.0,0.0),U∞=(0.0,0.0)))
    @test isapprox(Px, π*Rc^2, rtol=1e-1)
    @test isapprox(Py, 0.0, atol=1e-1)
    Px, Py = computeimpulse(model,parameters=ModelParameters(Ub=(0.0,1.0),U∞=(0.0,0.0)))
    @test isapprox(Px, 0.0, atol=1e-1)
    @test isapprox(Py, π*Rc^2, rtol=1e-1)

    model = VortexModel(g,bodies=ellipse)
    Px, Py = computeimpulse(model,parameters=ModelParameters(Ub=(1.0,0.0),U∞=(0.0,0.0)))
    @test isapprox(Px, π*b^2, rtol=1e-1)
    @test isapprox(Py, 0.0, atol=1e-1)
    Px, Py = computeimpulse(model,parameters=ModelParameters(Ub=(0.0,1.0),U∞=(0.0,0.0)))
    @test isapprox(Px, 0.0, atol=1e-1)
    @test isapprox(Py, π*a^2, rtol=1e-1)
end
