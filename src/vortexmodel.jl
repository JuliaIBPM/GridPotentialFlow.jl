import CartesianGrids: curl!
import LinearAlgebra: Diagonal

struct VortexModel{Nb,Ne,TU,TF}
    g::PhysicalGrid
    vortices::VortexList
    bodies::BodyList
    edges::Vector{Int}
    nodedata::TU
    bodydata::TF
    system::PotentialFlowSystem
end

function VortexModel(g::PhysicalGrid; bodies::Union{Body,Vector{<:Body},BodyList}=BodyList(), vortices::Union{Vortex,Vector{<:Vortex},VortexList}=VortexList(), edges::Vector{<:Integer}=Int[])

    # Ensure that bodies and vortices are of type BodyList and VortexList
    if isa(bodies,Body)
        bodies = BodyList([bodies])
    elseif isa(bodies,Vector{<:Body})
        bodies = BodyList(bodies)
    end
    if isa(vortices,Vortex)
        vortices = VortexList([vortices])
    elseif isa(vortices,Vector{<:Vortex})
        vortices = VortexList(vortices)
    end

    nodedata = Nodes(Dual,size(g))
    bodydata = ScalarData(length(collect(bodies)[1]))

    L = plan_laplacian(size(nodedata),with_inverse=true)

    if isempty(bodies) # Unregularized potential flow system without bodies
        system = PotentialFlowSystem(L)
    elseif isempty(edges) # Unregularized potential flow system with bodies
        regop = Regularize(VectorData(collect(bodies)),cellsize(g),I0=origin(g),issymmetric=true)
        Rmat,Emat = RegularizationMatrix(regop,bodydata,nodedata);
        S = SaddleSystem(L,Emat,Rmat,SaddleVector(nodedata,bodydata))
        system = PotentialFlowSystem(S)
    else  # Regularized potential flow system
        regop = Regularize(VectorData(collect(bodies)),cellsize(g),I0=origin(g),issymmetric=true)
        Rmat,Emat = RegularizationMatrix(regop,bodydata,nodedata);
        S = SaddleSystem(L,Emat,Rmat,SaddleVector(nodedata,bodydata))

        nodedata .= 0
        bodydata .= 1
        f₀ = constraint(S\SaddleVector(nodedata,bodydata));

        Df₀ = Diagonal(f₀);
        R̃mat = deepcopy(Rmat);
        R̃mat.M .= R̃mat.M*Df₀;
        S̃ = SaddleSystem(L,Emat,R̃mat,SaddleVector(nodedata,bodydata))

        e_kvec = [BodyUnitVector(bodies[1],k) for k in edges]
        d_kvec = typeof(nodedata)[]

        system = PotentialFlowSystem(S̃,f₀,e_kvec,d_kvec)
    end

    return VortexModel{length(bodies),length(edges),typeof(nodedata),typeof(bodydata)}(g, vortices, bodies, edges, nodedata, bodydata, system)
end

function issteady(vortexmodel::VortexModel)
    if isempty(vortexmodel.system.d_kvec)
        println("d_kvec not set in system. Model is treated as steady.")
        return true
    else
        return false
    end
end

function setvortices!(vortexmodel::VortexModel{Nb,Ne,TU,TF},vortices::Union{Vortex,Vector{<:Vortex},VortexList}=VortexList()) where {Nb,Ne,TU,TF}
    if isa(vortices,Vortex)
        vortices = VortexList([vortices])
    elseif isa(vortices,Vector{<:Vortex})
        vortices = VortexList(vortices)
    end

    vortexmodel.vortices = vortices
end

function setvortexpositions!(vortexmodel::VortexModel{Nb,Ne,TU,TF},X_vortices::VectorData{Nv}) where {Nb,Ne,TU,TF,Nv}

    @assert Nv == length(vortexmodel.vortices)

    setpositions!(vortexmodel.vortices,X_vortices.u,X_vortices.v)
end

function getvortexpositions(vortexmodel::VortexModel{Nb,Ne,TU,TF}) where {Nb,Ne,TU,TF}

    return getpositions(vortexmodel.vortices)
end

function step!(vortexmodel::VortexModel{Nb,Ne,TU,TF}) where {Nb,Ne,TU,TF} end

function computevortexvelocities(vortexmodel::VortexModel{Nb,Ne,TU,TF},X_vortices::VectorData{Nv};U∞=(0.0,0.0)) where {Nb,Ne,TU,TF,Nv}

    @assert Nv == length(vortexmodel.vortices)

    setvortexpositions!(vortexmodel,X_vortices)
    Ẋ_vortices = computevortexvelocities(vortexmodel,U∞=U∞)

    return Ẋ_vortices
end

function computevortexvelocities(vortexmodel::VortexModel{Nb,Ne,TU,TF}; U∞=(0.0,0.0)) where {Nb,Ne,TU,TF}

    Rmat,Emat = computeregularizationmatrix(vortexmodel)
    w = computew(vortexmodel,Rmat)
    ψ = computeψ(vortexmodel,w,U∞=U∞)
    Ẋ_vortices = computevortexvelocities(vortexmodel,ψ,Emat)

    return Ẋ_vortices
end

# function computevortexvelocities(vortexmodel::VortexModel{Nb,Ne,TU,TF},ψ::TU,Emat) where {Nb,Ne,TU,TF}
#
#     @unpack g, vortices = vortexmodel
#
#     q = NodePair(Dual,Dual,size(g))
#     curl!(q,ψ)
#
#     Ẋ = VectorData(length(collect(vortices)[1]));
#
#     Ẋ.u .= Emat*q.u
#     Ẋ.v .= Emat*q.v
#
#     return Ẋ
#
# end

function computevortexvelocities(vortexmodel::VortexModel{Nb,Ne,TU,TF},ψ::TU,Emat) where {Nb,Ne,TU,TF}

    @unpack g, vortices = vortexmodel

    Ẋ_vortices = VectorData(length(collect(vortices)[1]));
    s = TU();

    q = curl(ψ)/cellsize(g)

    grid_interpolate!(s,q.u);
    Ẋ_vortices.u .= Emat*s
    grid_interpolate!(s,q.v);
    Ẋ_vortices.v .= Emat*s

    return Ẋ_vortices
end

function computeregularizationmatrix(vortexmodel::VortexModel{Nb,Ne,TU,TF}) where {Nb,Ne,TU,TF}

    @unpack g, vortices, bodies, nodedata, bodydata = vortexmodel

    regop = Regularize(getpositions(vortices),cellsize(g),I0=origin(g),issymmetric=true)
    Rmat,Emat = RegularizationMatrix(regop,getstrengths(vortices), nodedata);

    return Rmat,Emat
end

function computew(vortexmodel::VortexModel{Nb,Ne,TU,TF}, Rmat::RegularizationMatrix{TU})::TU where {Nb,Ne,TU,TF}

    @unpack vortices = vortexmodel

    if isempty(vortices)
        return TU()
    end

    Γ = getstrengths(vortices)
    w = Rmat*Γ

    return w
end

function computew(vortexmodel::VortexModel{Nb,Ne,TU,TF})::TU where {Nb,Ne,TU,TF}

    @unpack g, vortices, bodies, nodedata, bodydata = vortexmodel

    if isempty(vortices)
        return TU()
    end

    Rmat,_ = computeregularizationmatrix(vortexmodel)
    w = computew(vortexmodel,Rmat)

    return w
end

function computeψ(vortexmodel::VortexModel{Nb,0,TU,TF}, w::TU; U∞=(0.0,0.0))::TU where {Nb,TU,TF}

    @unpack g, bodies, nodedata, bodydata, system = vortexmodel

    xg,yg = coordinates(nodedata,g)
    bodydata .= -U∞[1]*(collect(bodies)[2]) .+ U∞[2]*(collect(bodies)[1]);
    rhs = PotentialFlowRHS(w,bodydata)
    sol = PotentialFlowSolution(nodedata,bodydata)
    ldiv!(sol,system,rhs)
    sol.ψ .+= U∞[1]*yg' .- U∞[2]*xg

    return sol.ψ
end

function computeψ(vortexmodel::VortexModel{Nb,Ne,TU,TF}, w::TU; U∞=(0.0,0.0), σ=SuctionParameter.(zeros(Ne)))::TU where {Nb,Ne,TU,TF}

    @unpack g, bodies, edges, nodedata, bodydata, system = vortexmodel

    @assert length(edges) == length(σ)
    rhs = PotentialFlowRHS(w,bodydata,σ)

    xg,yg = coordinates(nodedata,g)
    bodydata .= -U∞[1]*(collect(bodies)[2]) .+ U∞[2]*(collect(bodies)[1]);

    # This line does the same as the next block of code, but creates new instances of TU and TF for the solution in systems.jl
    # sol = system\rhs

    if issteady(vortexmodel)
        sol = PotentialFlowSolution(nodedata,bodydata,zeros(Nb))
    else
        sol = PotentialFlowSolution(nodedata,bodydata,zeros(Nb),zeros(Ne))
    end
    ldiv!(sol,system,rhs)

    sol.ψ .+= U∞[1]*yg' .- U∞[2]*xg

    return sol.ψ
end

function computeψ(vortexmodel::VortexModel{Nb,Ne,TU,TF}; U∞=(0.0,0.0))::TU where {Nb,Ne,TU,TF}

    w = computew(vortexmodel)
    ψ = computeψ(vortexmodel, w, U∞=U∞)

    return ψ
end

# function computeψ(vortexmodel::VortexModel{Nb,Ne,TU,TF}; U∞=(0.0,0.0))::TU  where {Nb,Ne,TU,TF}
#
#     @unpack g, bodies, nodedata, bodydata, system = vortexmodel
#
#     xg,yg = coordinates(nodedata,g)
#     bodydata .= -U∞[1]*(collect(bodies)[2]) .+ U∞[2]*(collect(bodies)[1]);
#     rhs = PotentialFlowRHS(w,bodydata)
#     sol = PotentialFlowSolution(nodedata,bodydata,zeros(Nb),)
#
#     return sol.ψ
# end

function curl!(nodepair::NodePair{Dual, Dual, NX, NY},
               s::Nodes{Dual,NX, NY}) where {NX, NY}

    view(nodepair.u,2:NX-1,2:NY-1) .= 0.5*(view(s,2:NX-1,3:NY) .- view(s,2:NX-1,1:NY-2))
    #@inbounds for y in 1:NY-1, x in 1:NX
    #    edges.u[x,y] = s[x,y+1] - s[x,y]
    #end

    view(nodepair.v,2:NX-1,2:NY-1) .= 0.5*(view(s,1:NX-2,2:NY-1) .- view(s,3:NX,2:NY-1))
    #@inbounds for y in 1:NY, x in 1:NX-1
    #    edges.v[x,y] = s[x,y] - s[x+1,y]
    #end
    nodepair
end
