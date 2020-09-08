import CartesianGrids: curl!, Laplacian
import LinearAlgebra: Diagonal

mutable struct VortexModel{Nb,Ne,TU,TF}
    g::PhysicalGrid
    vortices::VortexList
    bodies::BodyList
    edges::Vector{Int}
    system::Union{PotentialFlowSystem,Laplacian}

    _nodedata::TU
    _bodydata::TF
    _w::TU
    _ψb::TF
    _Rmat::Union{Nothing,RegularizationMatrix{TU}}
    _Emat::Union{Nothing,InterpolationMatrix{TU}}
end

function VortexModel(g::PhysicalGrid; bodies::Union{Body,Vector{<:Body},BodyList}=BodyList(),  vortices::Union{Vortex,Vector{<:Vortex},VortexList}=VortexList(), edges::Vector{<:Integer}=Int[])

    # Ensure that bodies and vortices are of type BodyList and VortexList
    if isa(bodies,Body)
        bodies = BodyList([bodies])
    elseif isa(bodies,Vector{<:Body})
        bodies = BodyList(bodies)
    end

    _nodedata = Nodes(Dual,size(g))
    _bodydata = ScalarData(length(collect(bodies)[1]))
    _w = Nodes(Dual,size(g))
    _ψb = ScalarData(length(collect(bodies)[1]))

    L = plan_laplacian(size(_nodedata),with_inverse=true)
    regop = Regularize(VectorData(collect(bodies)), cellsize(g), I0=origin(g), issymmetric=true, ddftype=CartesianGrids.Yang3)
    Rmat,_ = RegularizationMatrix(regop,_bodydata,_nodedata);
    Emat = InterpolationMatrix(regop,_nodedata,_bodydata);
    S = SaddleSystem(L,Emat,Rmat,SaddleVector(_nodedata,_bodydata))
    if isempty(edges) # Unregularized potential flow system with bodies
        system = PotentialFlowSystem(S)
    else # Regularized potential flow system
        _nodedata .= 0
        _bodydata .= 1
        f₀ = constraint(S\SaddleVector(_nodedata,_bodydata));

        Df₀ = Diagonal(f₀);
        R̃mat = deepcopy(Rmat);
        R̃mat.M .= R̃mat.M*Df₀;
        S̃ = SaddleSystem(L,Emat,R̃mat,SaddleVector(_nodedata,_bodydata))

        e_kvec = [BodyUnitVector(bodies[1],k) for k in edges]
        d_kvec = typeof(_nodedata)[]

        system = PotentialFlowSystem(S̃,f₀,e_kvec,d_kvec)
    end

    vortexmodel =  VortexModel{length(bodies),length(edges),typeof(_nodedata),typeof(_bodydata)}(g, VortexList(), bodies, edges, system, _nodedata, _bodydata, _w, _ψb, nothing, nothing)

    setvortices!(vortexmodel, vortices)

    return vortexmodel
end

function issteady(vortexmodel::VortexModel)
    if isempty(vortexmodel.vortices)
        return true
    else
        return false
    end
end

function isregulated(vortexmodel::VortexModel{Nb,Ne}) where {Nb,Ne}
    if Ne > 0
        return true
    else
        return false
    end
end

function setvortices!(vortexmodel::VortexModel{Nb,Ne,TU,TF}, vortices::Union{Vortex,Vector{<:Vortex},VortexList}=VortexList(); computeregularization=true) where {Nb,Ne,TU,TF}

    vortices = VortexList(deepcopy(vortices))

    vortexmodel.vortices = vortices

    if computeregularization updatevorticesregularization(vortexmodel) end
end

function pushvortices!(vortexmodel::VortexModel{Nb,Ne,TU,TF}, vortices...; computeregularization=true) where {Nb,Ne,TU,TF}

    push!(vortexmodel.vortices.list,vortices...)

    if computeregularization updatevorticesregularization(vortexmodel) end
end

function setvortexpositions!(vortexmodel::VortexModel{Nb,Ne,TU,TF}, X_vortices::VectorData{Nv}; computeregularization=true) where {Nb,Ne,TU,TF,Nv}

    @assert Nv == length(vortexmodel.vortices)

    setpositions!(vortexmodel.vortices,X_vortices.u,X_vortices.v)

    if computeregularization updatevorticesregularization(vortexmodel) end
end

function getvortexpositions(vortexmodel::VortexModel{Nb,Ne,TU,TF}) where {Nb,Ne,TU,TF}

    return getpositions(vortexmodel.vortices)
end

function step!(vortexmodel::VortexModel{Nb,Ne,TU,TF}) where {Nb,Ne,TU,TF} end

function computevortexvelocities(vortexmodel::VortexModel{Nb,Ne,TU,TF},X_vortices::VectorData{Nv}; kwargs...) where {Nb,Ne,TU,TF,Nv}

    @assert Nv == length(vortexmodel.vortices)

    setvortexpositions!(vortexmodel,X_vortices)
    Ẋ_vortices = computevortexvelocities(vortexmodel; kwargs...)

    return Ẋ_vortices
end

function computevortexvelocities(vortexmodel::VortexModel{Nb,Ne,TU,TF}; kwargs...) where {Nb,Ne,TU,TF}

    computew!(vortexmodel._w,vortexmodel)
    sol = solvesystem!(vortexmodel._nodedata, vortexmodel._bodydata, vortexmodel, vortexmodel._w; kwargs...)

    for k in 1:Ne
        vortexmodel.vortices[end-Ne+k].Γ= sol.δΓ_kvec[k]/cellsize(vortexmodel.g)^2
    end

    Ẋ_vortices = computevortexvelocities(vortexmodel,sol.ψ)

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

function computevortexvelocities(vortexmodel::VortexModel{Nb,Ne,TU,TF},ψ::TU) where {Nb,Ne,TU,TF}

    @unpack g, vortices, _Emat = vortexmodel

    Ẋ_vortices = VectorData(length(collect(vortices)[1]));
    s = TU();

    q = curl(ψ)

    grid_interpolate!(s,q.u);
    Ẋ_vortices.u .= _Emat*s/cellsize(g)
    grid_interpolate!(s,q.v);
    Ẋ_vortices.v .= _Emat*s/cellsize(g)

    return Ẋ_vortices
end

function computeregularizationmatrix(g::PhysicalGrid,X::VectorData{N},f::ScalarData{N},s::Nodes) where {N}

    regop = Regularize(X, cellsize(g), I0=origin(g), ddftype=CartesianGrids.Yang3)
    Rmat = RegularizationMatrix(regop, f, s);
    Emat = InterpolationMatrix(regop, s, f);

    return Rmat, Emat
end

function updatevorticesregularization(vortexmodel::VortexModel{Nb,Ne,TU,TF}) where {Nb,Ne,TU,TF}

    @unpack g, vortices, _nodedata, system = vortexmodel

    vortexmodel._Rmat, vortexmodel._Emat = computeregularizationmatrix(vortexmodel.g,getpositions(vortices),getstrengths(vortices),vortexmodel._nodedata)

    if isregulated(vortexmodel) && !isempty(vortices)
        d_kvec = computed_kvec(vortexmodel,collect(length(vortices)-(Ne-1):length(vortices)))
        setd_kvec!(system, d_kvec)
    end
end

function computed_kvec(vortexmodel::VortexModel{Nb,Ne,TU,TF},indices::Vector{Int}) where {Nb,Ne,TU,TF}
    Γ = getstrengths(vortexmodel.vortices)
    d_kvec = TU[]
    for k in indices
        Γ .= 0
        Γ[k] = 1
        push!(d_kvec,vortexmodel._Rmat*Γ)
    end
    return d_kvec
end

function computew(vortexmodel::VortexModel{Nb,Ne,TU,TF})::TU where {Nb,Ne,TU,TF}

    @unpack g, vortices, _Rmat = vortexmodel

    w = TU()
    computew!(w,vortexmodel)

    return w
end

function computew!(w,vortexmodel::VortexModel{Nb,Ne,TU,TF})::TU where {Nb,Ne,TU,TF}

    @unpack g, vortices, _Rmat = vortexmodel

    if isempty(vortices)
        w .= 0.0
        return w
    end

    Γ = getstrengths(vortices)
    Γ[end-Ne+1:end] .= 0
    w .= _Rmat*Γ # _Rmat takes physical cell size into account

    return w
end

function solvesystem!(sol, vortexmodel::VortexModel{Nb,Ne,TU,TF}, wphysical::TU; Ub=(0.0,0.0), U∞=(0.0,0.0), Γb=nothing, σ=SuctionParameter.(zeros(Ne))) where {Nb,Ne,TU,TF}

    @unpack g, vortices, bodies, edges, system, _nodedata, _bodydata, _w, _ψb = vortexmodel

    # Scale w to grid with unit cell size
    _w .= wphysical*cellsize(g)^2

    _ψb .= (Ub[1]-U∞[1])*(collect(bodies)[2]) .+ (-Ub[2]+U∞[2])*(collect(bodies)[1]);

    if !isregulated(vortexmodel)
        rhs = PotentialFlowRHS(_w,_ψb,Γb)
    else
        rhs = PotentialFlowRHS(_w,_ψb,σ)
    end

    ldiv!(sol,system,rhs)

    xg,yg = coordinates(_nodedata,g)
    sol.ψ .+= U∞[1]*yg' .- U∞[2]*xg

    return sol
end

function solvesystem!(ψ::TU, f::TF, vortexmodel::VortexModel{Nb,Ne,TU,TF}, wphysical::TU; kwargs...) where {Nb,Ne,TU,TF}

    if !isregulated(vortexmodel)
        sol = PotentialFlowSolution(ψ,f)
    elseif issteady(vortexmodel)
        sol = PotentialFlowSolution(ψ,f,zeros(Nb))
    else
        sol = PotentialFlowSolution(ψ,f,zeros(Nb),zeros(Ne))
    end

    solvesystem!(sol,vortexmodel,wphysical;kwargs...)

    return sol
end

function solvesystem(vortexmodel::VortexModel{Nb,Ne,TU,TF}, wphysical::TU; kwargs...) where {Nb,Ne,TU,TF}

    solvesystem!(TU(),TF(),vortexmodel,wphysical;kwargs...)

    return sol
end

function computeψ(vortexmodel::VortexModel{Nb,Ne,TU,TF}; kwargs...)::TU where {Nb,Ne,TU,TF}

    ψ = TU()

    computew!(vortexmodel._w, vortexmodel)
    solvesystem!(ψ, vortexmodel._bodydata, vortexmodel, vortexmodel._w; kwargs...)

    return ψ
end

# function computeψ(vortexmodel::VortexModel{Nb,Ne,TU,TF}; U∞=(0.0,0.0))::TU  where {Nb,Ne,TU,TF}
#
#     @unpack g, bodies, system, _nodedata, _bodydata = vortexmodel
#
#     xg,yg = coordinates(_nodedata,g)
#     _bodydata .= -U∞[1]*(collect(bodies)[2]) .+ U∞[2]*(collect(bodies)[1]);
#     rhs = PotentialFlowRHS(w,_bodydata)
#     sol = PotentialFlowSolution(_nodedata,_bodydata,zeros(Nb),)
#
#     return sol.ψ
# end

function computeimpulse(vortexmodel::VortexModel{Nb,Ne,TU,TF}, w::TU, f::TF; Ub=(0.0,0.0), U∞=(0.0,0.0), kwargs...) where {Nb,Ne,TU,TF}

    @unpack g, vortices, bodies = vortexmodel

    xg, yg = coordinates(w,g)
    Δx = cellsize(g)

    # Formula 6.16
    volumeintegral_x = Δx^2*sum(w.*yg)
    volumeintegral_y = Δx^2*sum(-w.*xg)
    if !isempty(bodies)
        ncrossUb = -(normalmid(bodies[1])[1]*Ub[2] - normalmid(bodies[1])[2]*Ub[1])
        surfaceintegral_x = _surfaceintegrate(bodies[1],(f./dlength(bodies[1]) + ncrossUb).*(bodies[1].y))
        surfaceintegral_y = _surfaceintegrate(bodies[1],-(f./dlength(bodies[1]) + ncrossUb).*(bodies[1].x))
    end

    P_x = volumeintegral_x + surfaceintegral_x - _calculatevolume(bodies[1])*U∞[1]
    P_y = volumeintegral_y + surfaceintegral_y - _calculatevolume(bodies[1])*U∞[2]

    return P_x, P_y
end

function computeimpulse(vortexmodel::VortexModel{Nb,Ne,TU,TF}; kwargs...) where {Nb,Ne,TU,TF}

    computew!(vortexmodel._w,vortexmodel)
    solvesystem!(vortexmodel._nodedata, vortexmodel._bodydata, vortexmodel, vortexmodel._w; kwargs...)

    P_x, P_y = computeimpulse(vortexmodel, vortexmodel._w, vortexmodel._bodydata; kwargs...)

    return P_x, P_y
end

# function curl!(nodepair::NodePair{Dual, Dual, NX, NY},
#                s::Nodes{Dual,NX, NY}) where {NX, NY}
#
#     view(nodepair.u,2:NX-1,2:NY-1) .= 0.5*(view(s,2:NX-1,3:NY) .- view(s,2:NX-1,1:NY-2))
#     #@inbounds for y in 1:NY-1, x in 1:NX
#     #    edges.u[x,y] = s[x,y+1] - s[x,y]
#     #end
#
#     view(nodepair.v,2:NX-1,2:NY-1) .= 0.5*(view(s,1:NX-2,2:NY-1) .- view(s,3:NX,2:NY-1))
#     #@inbounds for y in 1:NY, x in 1:NX-1
#     #    edges.v[x,y] = s[x,y] - s[x+1,y]
#     #end
#     nodepair
# end

# Only works for closed bodies for now
function _surfaceintegrate(body::Body{N,C},integrand::Array{Float64,1}) where {N,C}
    func = Array{Float64,1}(undef, N+1)
    func[1:end-1] .= integrand
    func[end] = integrand[1]
    s = sum(dlength(body).*(func[1:end-1] + func[2:end]))/2
    return s
end

function _calculatevolume(body::Body{N,RigidBodyTools.ClosedBody}) where {N}

    ip1(i) = 1+mod(i,N)
    V = sum([body.x[i]*body.y[ip1(i)] - body.y[i]*body.x[ip1(i)] for i = 1:N])/2
    return V
end
