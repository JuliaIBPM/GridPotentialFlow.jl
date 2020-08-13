import CartesianGrids: curl!

struct VortexModel{Nb,TU,TF}
    g::PhysicalGrid
    vortices::VortexList
    bodies::BodyList
    nodedata::TU
    bodydata::TF
    system::Union{CartesianGrids.Laplacian,PotentialFlowSystem}
end

function VortexModel(g::PhysicalGrid; bodies::Union{Body,Vector{<:Body},BodyList}=BodyList(), vortices::Union{Vortex,Vector{<:Vortex},VortexList}=VortexList(), edges::Vector{Float64}=Float64[])

    # Assert that bodies and vortices are of type BodyList and VortexList
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

    if isempty(bodies)
        system = L
    elseif isempty(edges)
        regop = Regularize(VectorData(collect(bodies)),cellsize(g),I0=origin(g),issymmetric=true)
        Rmat,Emat = RegularizationMatrix(regop,bodydata,nodedata);
        S = SaddleSystem(L,Emat,Rmat,SaddleVector(nodedata,bodydata))
        system = PotentialFlowSystem(S)
    else
        regop = Regularize(VectorData(collect(bodies)),cellsize(g),I0=origin(g),issymmetric=true)
        Rmat,Emat = RegularizationMatrix(regop,bodydata,nodedata);
        S = SaddleSystem(L,Emat,Rmat,SaddleVector(nodedata,bodydata))

        f₀ = constraint(S\SaddleVector(nodedata,bodydata));
        e_kvec = [BodyUnitVector(plate,k) for k in edges]
        d_kvec = [deepcopy(nodedata) for k in edges]

        system = PotentialFlowSystem(S,f₀,e_kvec,d_kvec)
    end

    return VortexModel{length(bodies),typeof(nodedata),typeof(bodydata)}(g, vortices, bodies, nodedata, bodydata, system)
end

function setvortices!(vortexmodel::VortexModel{Nb,TU,TF},vortices::Union{Vortex,Vector{<:Vortex},VortexList}=VortexList()) where {Nb,TU,TF}

    if isa(vortices,Vortex)
        vortices   = VortexList([vortices])
    elseif isa(vortices,Vector{<:Vortex})
        vortices = VortexList(vortices)
    end

    vortexmodel.vortices = vortices

end

function setvortexpositions!(vortexmodel::VortexModel{Nb,TU,TF},X_vortices::VectorData{Nv}) where {Nb,TU,TF,Nv}

    @assert Nv == length(vortexmodel.vortices)

    setpositions!(vortexmodel.vortices,X_vortices.u,X_vortices.v)

end

function getvortexpositions(vortexmodel::VortexModel{Nb,TU,TF}) where {Nb,TU,TF}

    getpositions(vortexmodel.vortices)

end

function step!(vortexmodel::VortexModel{Nb,TU,TF}) where {Nb,TU,TF} end

function computevelocity(vortexmodel::VortexModel{Nb,TU,TF},X_vortices::VectorData{Nv};U∞=(0.0,0.0)) where {Nb,TU,TF,Nv}

    @assert Nv == length(vortexmodel.vortices)

    setvortexpositions!(vortexmodel,X_vortices)
    Ẋ_vortices = computevelocity(vortexmodel,U∞=U∞)

    return Ẋ_vortices

end

function computevelocity(vortexmodel::VortexModel{Nb,TU,TF}; U∞=(0.0,0.0)) where {Nb,TU,TF}

    Rmat,Emat = computeregularizationmatrix(vortexmodel)
    w = computew(vortexmodel,Rmat)
    ψ = computeψ(vortexmodel,w,U∞=U∞)
    Ẋ_vortices = computevelocity(vortexmodel,ψ,Emat)

    return Ẋ_vortices
end

# function computevelocity(vortexmodel::VortexModel{Nb,TU,TF},ψ::TU,Emat) where {Nb,TU,TF}
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

function computevelocity(vortexmodel::VortexModel{Nb,TU,TF},ψ::TU,Emat) where {Nb,TU,TF}

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

function computeregularizationmatrix(vortexmodel::VortexModel{Nb,TU,TF}) where {Nb,TU,TF}

    @unpack g, vortices, bodies, nodedata, bodydata = vortexmodel

    regop = Regularize(getpositions(vortices),cellsize(g),I0=origin(g),issymmetric=true)
    Rmat,Emat = RegularizationMatrix(regop,getstrengths(vortices), nodedata);

    return Rmat,Emat
end

function computew(vortexmodel::VortexModel{Nb,TU,TF}, Rmat::RegularizationMatrix{TU})::TU where {Nb,TU,TF}

    @unpack vortices = vortexmodel

    if isempty(vortices)
        return TU()
    end

    Γ = getstrengths(vortices)
    w = Rmat*Γ

    return w
end

function computew(vortexmodel::VortexModel{Nb,TU,TF})::TU where {Nb,TU,TF}

    @unpack g, vortices, bodies, nodedata, bodydata = vortexmodel

    if isempty(vortices)
        return TU()
    end

    Rmat,_ = computeregularizationmatrix(vortexmodel)
    w = computew(vortexmodel,Rmat)

    return w
end

function computeψ(vortexmodel::VortexModel{0,TU,TF}, w::TU; U∞=(0.0,0.0))::TU where {Nb,TU,TF}

    @unpack g, bodies, nodedata, bodydata, system = vortexmodel

    xg,yg = coordinates(nodedata,g)
    ψ = TU()
    ldiv!(ψ,system,w)
    ψ .+= U∞[1]*yg' .- U∞[2]*xg

    return ψ
end

function computeψ(vortexmodel::VortexModel{Nb,TU,TF}, w::TU; U∞=(0.0,0.0))::TU where {Nb,TU,TF}

    @unpack g, bodies, nodedata, bodydata, system = vortexmodel

    xg,yg = coordinates(nodedata,g)
    bodydata .= -U∞[1]*(collect(bodies)[2]) .+ U∞[2]*(collect(bodies)[1]);
    rhs = PotentialFlowRHS(w,bodydata)
    sol = PotentialFlowSolution(nodedata,bodydata)
    ldiv!(sol,system,rhs)
    sol.ψ .+= U∞[1]*yg' .- U∞[2]*xg

    return sol.ψ
end

function computeψ(vortexmodel::VortexModel{Nb,TU,TF}; U∞=(0.0,0.0))::TU where {Nb,TU,TF}

    w = computew(vortexmodel)
    ψ = computeψ(vortexmodel, w, U∞=U∞)

    return ψ
end

# function computeψ(vortexmodel::VortexModel{Nb,TU,TF};U∞,edges,vortices)::TU  where {Nb,TU,TF} end

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
