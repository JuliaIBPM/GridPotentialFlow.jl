import CartesianGrids: curl!, Laplacian
import LinearAlgebra: Diagonal, norm

export VortexModel, computeψ, computew, computew!, computevortexvelocities, computeregularizationmatrix, getstrengths, getpositions, setvortexpositions!, getvortexpositions, setvortices!, pushvortices!, computeimpulse, computeaddedmassmatrix, solvesystem, solvesystem!

# TODO: check if _computeψboundaryconditions needs to be faster
# TODO: impulse case when U∞ is specified instead of Ub
# TODO: create PotentialFlowBody type that encapsulates Ub, Γb, and any regularized edges
# TODO: try to remove _d_kvec from VortexModel
# TODO: check memory allocation inverse laplacian in ConstrainedSystems
# TODO: match names of Ne and Nk

mutable struct VortexModel{Nb,Ne,isshedding}
    g::PhysicalGrid
    bodies::BodyList
    vortices::VortexList
    edges::Vector{Int}
    system::PotentialFlowSystem

    _nodedata::Nodes{Dual}
    _edgedata::Edges{Primal}
    _bodydata::ScalarData
    _bodyvectordata::VectorData
    _d_kvec::Vector{Nodes{Dual}}
    _ψ::Nodes{Dual}
    _f::ScalarData
    _w::Nodes{Dual}
    _ψb::ScalarData
end

function VortexModel(g::PhysicalGrid; bodies::Union{Body,Vector{<:Body},BodyList}=BodyList(), vortices::Union{Vortex,Vector{<:Vortex},VortexList}=VortexList(), edges::Vector{<:Integer}=Int[])

    # Ensure that bodies are of type BodyList
    # TODO: check other ways for code robustness
    if bodies isa Body
        bodies = BodyList([bodies])
    elseif bodies isa Vector{<:Body}
        bodies = BodyList(bodies)
    end

    Nb = length(bodies)
    Ne = length(edges)

    # Initialize data structures for internal use
    _nodedata = Nodes(Dual,size(g))
    _edgedata = Edges(Primal,size(g))
    _bodydata = ScalarData(length(collect(bodies)[1]))
    _bodyvectordata = VectorData(length(collect(bodies)[1]))
    _d_kvec = typeof(_nodedata)[]
    _ψ = Nodes(Dual,size(g))
    _f = ScalarData(length(collect(bodies)[1]))
    _w = Nodes(Dual,size(g))
    _ψb = ScalarData(length(collect(bodies)[1]))

    _TF_ones = zeros(length(collect(bodies)[1]),Nb)
    for i in 1:Nb
        _TF_ones[getrange(bodies,i),i] .= 1;
    end

    if Ne == 0 # Without Kutta condition
        L = plan_laplacian(size(_nodedata),with_inverse=true)
        _Rbmat, _Ebmat = computeregularizationmatrix(g,VectorData(collect(bodies)), _bodydata, _nodedata)
        S = SaddleSystem(L,_Ebmat,_Rbmat,SaddleVector(_nodedata,_bodydata))
        system = PotentialFlowSystem(S,_TF_ones)
    else # With Kutta condition
        _d_kvec = [Nodes(Dual,size(g)) for i=1:Ne]
        L = plan_laplacian(size(_nodedata),with_inverse=true)
        _Rbmat, _Ebmat = computeregularizationmatrix(g,VectorData(collect(bodies)), _bodydata, _nodedata)
        S = SaddleSystem(L,_Ebmat,_Rbmat,SaddleVector(_nodedata,_bodydata))
        _nodedata .= 0
        _bodydata .= 1
        f₀ = constraint(S\SaddleVector(_nodedata,_bodydata));
        R̃bmat = deepcopy(_Rbmat);
        R̃bmat.M .= R̃bmat.M*Diagonal(f₀);
        S̃ = SaddleSystem(L,_Ebmat,R̃bmat,SaddleVector(_nodedata,_bodydata))
        e_kvec = [BodyUnitVector(bodies[1],k) for k in edges]
        system = PotentialFlowSystem(S̃,f₀,e_kvec,_d_kvec)
    end

    if !isempty(vortices) && Ne > 0
        isshedding = true
    else
        isshedding = false
    end

    vortexmodel =  VortexModel{Nb,Ne,isshedding}(g, bodies, VortexList(deepcopy(vortices)), edges, system, _nodedata, _edgedata, _bodydata, _bodyvectordata, _d_kvec, _ψ, _f, _w, _ψb)

    return vortexmodel
end

function setvortices!(vortexmodel::VortexModel{Nb,Ne}, vortices::Union{Vortex,Vector{<:Vortex},VortexList}=VortexList()) where {Nb,Ne}

    vortices = VortexList(deepcopy(vortices))

    vortexmodel.vortices = vortices
end

function pushvortices!(vortexmodel::VortexModel{Nb,Ne}, vortices...) where {Nb,Ne}

    push!(vortexmodel.vortices.list,vortices...)
end

function setvortexpositions!(vortexmodel::VortexModel{Nb,Ne}, X_vortices::VectorData{Nv}) where {Nb,Ne,Nv}

    @assert Nv == length(vortexmodel.vortices)

    setpositions!(vortexmodel.vortices,X_vortices.u,X_vortices.v)
end

function getvortexpositions(vortexmodel::VortexModel{Nb,Ne}) where {Nb,Ne}

    return getpositions(vortexmodel.vortices)
end

function _updatesystemd_kvec!(vortexmodel::VortexModel{Nb,Ne}, indices::Vector{Int}) where {Nb,Ne}

    @unpack g, vortices, system, _nodedata, _d_kvec = vortexmodel

    H = Regularize(getpositions(vortices), cellsize(g), I0=origin(g), ddftype=CartesianGrids.M4prime, issymmetric=true)

    Γ = getstrengths(vortices)
    for i in 1:Ne
        Γ .= 0
        Γ[indices[i]] = 1
        H(_nodedata,Γ)
        _d_kvec[i] .= _nodedata
    end
    setd_kvec!(system,_d_kvec)
    return _d_kvec
end

"""
    computevortexvelocities(vortexmodel::VortexModel, ψ::Nodes{Dual})

Returns the flow velocity as `VectorData` at the locations of the vortices stored in `vortexmodel` associated with the discrete vector potential field `ψ`.
"""
function computevortexvelocities(vortexmodel::VortexModel{Nb,Ne}, ψ::Nodes{Dual}) where {Nb,Ne}

    @unpack g, vortices, _nodedata, _edgedata, _bodydata = vortexmodel

    H = Regularize(getpositions(vortices), cellsize(g), I0=origin(g), ddftype=CartesianGrids.M4prime, issymmetric=true)

    Ẋ_vortices = VectorData(length(vortices))
    # Velocity is the curl of the vector potential
    # The discrete curl operator requires dividing by the cellsize to account for the grid spacing
    curl!(_edgedata,ψ)
    _edgedata ./= cellsize(g)

    # For consistent interpolation, first interpolate the velocity to the nodes and use H to interpolate from the nodes to the vortices
    grid_interpolate!(_nodedata,_edgedata.u);
    H(Ẋ_vortices.u,_nodedata)
    grid_interpolate!(_nodedata,_edgedata.v);
    H(Ẋ_vortices.v,_nodedata)

    return Ẋ_vortices
end

"""
    computevortexvelocities(vortexmodel::VortexModel, w::Nodes{Dual}; kwargs...)

Returns the flow velocity as `VectorData` at the locations of the vortices stored in `vortexmodel` due to the vorticity field `w` and accounting for bodies in 'vortexmodel' and conditions in 'kwargs'.
"""
function computevortexvelocities(vortexmodel::VortexModel{Nb,Ne}; kwargs...) where {Nb,Ne}

    @unpack g, vortices, _nodedata, _bodydata = vortexmodel

    # The strengths of the Ne last vortices will be calculated in solvesystem and should be set to zero before computing the vorticity field such that they are not included in w
    for k in 1:Ne
        vortexmodel.vortices[end-Ne+k].Γ = 0.0
    end

    computew!(vortexmodel._w,vortexmodel)

    sol = solvesystem(vortexmodel, vortexmodel._w; kwargs...)

    for k in 1:Ne
        vortices[end-Ne+k].Γ = sol.δΓ_kvec[k]
    end

    Ẋ_vortices = computevortexvelocities(vortexmodel,sol.ψ)

    return Ẋ_vortices
end

function computevortexvelocities(vortexmodel::VortexModel{Nb,Ne}, X_vortices::VectorData{Nv}; kwargs...) where {Nb,Ne,Nv}

    @assert Nv == length(vortexmodel.vortices)

    setvortexpositions!(vortexmodel,X_vortices)
    Ẋ_vortices = computevortexvelocities(vortexmodel; kwargs...)

    return Ẋ_vortices
end

"""
    computew!(w::Nodes,vortexmodel::VortexModel)

Computes the vorticity field `w` associated with the vortices stored in `vortexmodel` on the physical grid.
"""
function computew!(wphysical::Nodes{Dual}, vortexmodel::VortexModel{Nb,Ne})::Nodes{Dual} where {Nb,Ne}

    @unpack g, vortices = vortexmodel

    H = Regularize(getpositions(vortices), cellsize(g), I0=origin(g), ddftype=CartesianGrids.M4prime, issymmetric=true)

    if isempty(vortices)
        wphysical .= 0.0
        return wphysical
    end

    Γ = getstrengths(vortices)
    H(wphysical,Γ)
    wphysical ./= cellsize(g)^2 # Divide by the Δx² to ensure that ∫wdA = ΣΓ

    return wphysical
end

"""
    computew(sol::PotentialFlowSolution, vortexmodel::VortexModel, wphysical::Nodes, parameters::ModelParameters)

"""
function computew(vortexmodel::VortexModel{Nb,Ne})::Nodes{Dual} where {Nb,Ne}

    @unpack g, vortices = vortexmodel

    w = Nodes(Dual,size(g))
    computew!(w,vortexmodel)

    return w
end

"""
    solvesystem!(sol::PotentialFlowSolution, vortexmodel::VortexModel, wphysical::Nodes, parameters::ModelParameters)

"""
function solvesystem!(sol::UnregularizedPotentialFlowSolution, vortexmodel::VortexModel{Nb,0}, wphysical::Nodes{Dual}; parameters=ModelParameters()) where {Nb}

    _computeψboundaryconditions!(vortexmodel._ψb, vortexmodel, parameters)
    vortexmodel._w .= wphysical
    Γb = deepcopy(parameters.Γb)

    if isnothing(Γb) && Nb == 1 # Circulation about body not specified, using Γb = -∫wdA
        Γb = -sum(vortexmodel._w)*cellsize(vortexmodel.g)^2
    elseif isnothing(Γb) && Nb > 1 # Circulation about bodies not specified, using Γb = zeros(Nb)
        Γb = zeros(Nb)
    end

    rhs = PotentialFlowRHS(vortexmodel._w,vortexmodel._ψb,Γ=Γb)
    _scaletoindexspace!(rhs,cellsize(vortexmodel.g))
    ldiv!(sol,vortexmodel.system,rhs)
    _scaletophysicalspace!(sol,cellsize(vortexmodel.g))

    addψ∞!(sol.ψ,vortexmodel,parameters) # Add the uniform flow to the approximation to the continuous stream function field

    return sol
end

function solvesystem!(sol::SteadyRegularizedPotentialFlowSolution, vortexmodel::VortexModel{Nb,Ne,false}, wphysical::Nodes{Dual}; parameters=ModelParameters()) where {Nb,Ne}

    _computeψboundaryconditions!(vortexmodel._ψb, vortexmodel, parameters)

    vortexmodel._w .= wphysical
    SP = isnothing(parameters.σ) ? SuctionParameter.(zeros(Ne)) : deepcopy(parameters.σ)

    rhs = PotentialFlowRHS(vortexmodel._w,vortexmodel._ψb,SP)
    _scaletoindexspace!(rhs,cellsize(vortexmodel.g))
    ldiv!(sol,vortexmodel.system,rhs)
    _scaletophysicalspace!(sol,cellsize(vortexmodel.g))

    addψ∞!(sol.ψ,vortexmodel,parameters) # Add the uniform flow to the approximation to the continuous stream function field

    return sol
end

function solvesystem!(sol::UnsteadyRegularizedPotentialFlowSolution, vortexmodel::VortexModel{Nb,Ne,true}, wphysical::Nodes{Dual}; parameters=ModelParameters()) where {Nb,Ne}

    _computeψboundaryconditions!(vortexmodel._ψb, vortexmodel, parameters)
    _updatesystemd_kvec!(vortexmodel,collect(length(vortexmodel.vortices)-(Ne-1):length(vortexmodel.vortices))) # Update the system.d_kvec to agree with the latest vortex positions

    vortexmodel._w .= wphysical
    SP = isnothing(parameters.σ) ? SuctionParameter.(zeros(Ne)) : deepcopy(parameters.σ)
    Γw = sum(vortexmodel._w)*cellsize(vortexmodel.g)^2 # Γw = ∫wdA

    rhs = PotentialFlowRHS(vortexmodel._w,vortexmodel._ψb,SP,Γw)
    _scaletoindexspace!(rhs,cellsize(vortexmodel.g))
    ldiv!(sol,vortexmodel.system,rhs)
    _scaletophysicalspace!(sol,cellsize(vortexmodel.g))

    addψ∞!(sol.ψ,vortexmodel,parameters) # Add the uniform flow to the approximation to the continuous stream function field

    return sol
end

"""
    solvesystem(vortexmodel::VortexModel, wphysical::Nodes; parameters::ModelParameters)

"""
function solvesystem(vortexmodel::VortexModel{Nb,0,false}, wphysical::Nodes{Dual}; kwargs...) where {Nb}

    sol = PotentialFlowSolution(vortexmodel._ψ, vortexmodel._f)
    solvesystem!(sol, vortexmodel, wphysical; kwargs...)

    return sol
end

function solvesystem(vortexmodel::VortexModel{Nb,Ne,false}, wphysical::Nodes{Dual}; kwargs...) where {Nb,Ne}

    sol = PotentialFlowSolution(vortexmodel._ψ, vortexmodel._f, zeros(Nb))
    solvesystem!(sol, vortexmodel, wphysical; kwargs...)

    return sol
end

function solvesystem(vortexmodel::VortexModel{Nb,Ne,true}, wphysical::Nodes{Dual}; kwargs...) where {Nb,Ne}

    sol = PotentialFlowSolution(vortexmodel._ψ, vortexmodel._f, zeros(Nb), zeros(Ne))
    solvesystem!(sol, vortexmodel, wphysical; kwargs...)

    return sol
end

function computeψ(vortexmodel::VortexModel; kwargs...)::Nodes{Dual}

    computew!(vortexmodel._w, vortexmodel)
    sol = solvesystem(vortexmodel, vortexmodel._w; kwargs...)

    return sol.ψ
end

function addψ∞!(ψ, vortexmodel::VortexModel, parameters)::Nodes{Dual}

    xg,yg = coordinates(vortexmodel._nodedata,vortexmodel.g)
    ψ .+= parameters.U∞[1].*yg' .- parameters.U∞[2].*xg

    return vortexmodel._nodedata
end

function _computeψboundaryconditions!(ψb::ScalarData, vortexmodel::VortexModel, parameters)

    _computebodypointsvelocity!(vortexmodel._bodyvectordata, parameters.Ub, vortexmodel.bodies) # Convert Ub into VectorData corresponding with the body points

    # The discrete streamfunction field is constrained to a prescribed streamfunction on the body that describes the body motion. The body presence in the uniform flow is taken into account by subtracting its value from the body motion (i.e. a body motion in the -U∞ direction) and adding the uniform flow at the end of the solvesystem routine.
    ψb .= -parameters.U∞[1]*(collect(vortexmodel.bodies)[2]) .+ parameters.U∞[2]*(collect(vortexmodel.bodies)[1]);
    ψb .+= vortexmodel._bodyvectordata.u .* collect(vortexmodel.bodies)[2] .- vortexmodel._bodyvectordata.v .* collect(vortexmodel.bodies)[1]

    return ψb
end

function _computeψboundaryconditions!(ψb::ScalarData, vortexmodel::VortexModel{0}, parameters)
    return
end

# frame of reference

# w and f are physical quantities
function computeimpulse(vortexmodel::VortexModel{Nb,Ne}, w::Nodes{Dual}, f::ScalarData, Ubvec) where {Nb,Ne}

    @unpack g, vortices, bodies, _bodydata, _bodyvectordata = vortexmodel

    xg, yg = coordinates(w,g)
    Δx = cellsize(g)

    # Formula 61 (see formula 6.16 in book)
    impulse = [Δx^2*sum(w.*yg'),Δx^2*sum(-w.*xg)]

    for i in 1:Nb
        impulse += computeimpulsesurfaceintegral(bodies[i], f[getrange(bodies,i)], Ubvec.u[getrange(bodies,i)], Ubvec.v[getrange(bodies,i)])
    end

    return impulse[1], impulse[2]
end

function computeimpulse(vortexmodel::VortexModel; parameters=ModelParameters()) where {Nb}

    @unpack bodies, _nodedata, _bodydata, _bodyvectordata, _w = vortexmodel

    computew!(_nodedata,vortexmodel)
    # Solve system, _bodydata will contain f
    sol = solvesystem(vortexmodel, _nodedata, parameters=parameters)

    # Convert Ub into VectorData corresponding with the body points
    _computebodypointsvelocity!(_bodyvectordata,parameters.Ub,bodies)

    P_x, P_y = computeimpulse(vortexmodel, _nodedata, sol.f, _bodyvectordata)

    return P_x, P_y
end

function computeimpulsesurfaceintegral(body::Body{N,RigidBodyTools.ClosedBody}, f, u, v) where {N}
    nx,ny = normalmid(body)
    Δs = dlengthmid(body)
    return computecrossproductsurfaceintegral(body,(f./Δs + nx.*v - ny.*u).*Δs)
end

function computeimpulsesurfaceintegral(body::Body{N,RigidBodyTools.OpenBody}, f, u, v) where {N}
    return computecrossproductsurfaceintegral(body,f)
end

function computeaddedmassmatrix(vortexmodel::VortexModel{Nb,Ne}) where {Nb,Ne}
    @unpack g, vortices, bodies, _bodyvectordata, _w = vortexmodel

    # For now only translational
    M = zeros(Nb*2,Nb*2)

    computew!(_w,vortexmodel)

    for movingbodyindex in 1:Nb
        for dir in 1:2
            Ub = fill((0.0,0.0),Nb)
            if dir == 1
                Ub[movingbodyindex] = (1.0,0.0)
            else
                Ub[movingbodyindex] = (0.0,1.0)
            end
            _computebodypointsvelocity!(_bodyvectordata,Ub,bodies)
            sol = solvesystem(vortexmodel, _w; parameters=ModelParameters(Ub=Ub))
            for i in 1:Nb
                M[(i-1)*2+1:(i-1)*2+2,(movingbodyindex-1)*2+dir] .= computeimpulsesurfaceintegral(bodies[i], sol.f[getrange(bodies,i)], _bodyvectordata.u[getrange(bodies,i)], _bodyvectordata.v[getrange(bodies,i)])
            end
        end
    end
    return M
end
