using GridPotentialFlow
using Test
using LinearAlgebra: eigen

function rectangulararray(rows,columns,spacing)
    centers = Tuple{Float64,Float64}[]
    for j in 1:columns
        for i in 1:rows
            xc = (i-1)*spacing - (columns-1)*spacing/2
            yc = (j-1)*spacing - (rows-1)*spacing/2
            push!(centers,(xc,yc))
        end
    end
    return centers
end

@testset "Added mass" begin
    Δx = 0.01
    xlim = (-2,2)
    ylim = (-2,2)
    g = PhysicalGrid(xlim,ylim,Δx)
    Δs = 1.4*cellsize(g)

    @testset "Ellipse" begin
        a = 0.5
        b = 2
        ellipse = Ellipse(a,b,Δs)
        prob = setup_problem(g,ellipse,phys_params=Dict())
        sys = construct_system(prob)
        M = addedmass(sys)
        @test isapprox(M[1,1], π*b^2, rtol=1e-1)
        @test isapprox(M[2,2], π*a^2, rtol=1e-1)
    end

    @testset "Cylinder array" begin
        nineRodArrayChen = [2.3637,2.2092,2.1007,2.0120,1.9350,1.8665,1.7494,1.6531,1.5732,1.5066,1.4508]
        GRratios = [0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.2,1.4,1.6]
        idx = length(GRratios)
        rows = 3
        columns = 3
        N = rows*columns
        R = 0.20
        𝒱 = π*R^2
        bodylist = BodyList([Circle(R,Δs) for i in 1:N])
        gap = GRratios[idx]*R
        spacing = 2*R+gap
        bodycenters = rectangulararray(rows,columns,spacing)
        tflist = RigidTransformList(RigidTransform.(bodycenters,0.0))
        tflist(bodylist)
        prob = PotentialFlowProblem(g,bodylist)
        sys = construct_system(prob)
        M = addedmass(sys)/𝒱
        eigM = eigen(M);
        max_eig_value_coef = maximum(real(eigM.values))
        max_self_added_mass_coef = maximum(M)
        λoverMratio = max_eig_value_coef/max_self_added_mass_coef
        @test isapprox(λoverMratio,nineRodArrayChen[idx],rtol=0.01)
    end

end
