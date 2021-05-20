import CartesianGrids: curl!, Laplacian
import LinearAlgebra: Diagonal, norm
import Base: show

export VortexModel, computeψ, computew, computew!, computevortexvelocities, _computeregularizationmatrix, getstrengths, getpositions, setvortexpositions!, getvortexpositions, setvortices!, pushvortices!, computeimpulse, computeaddedmassmatrix, solvesystem, solvesystem!

# TODO: check if _computeψboundaryconditions needs to be faster
# TODO: impulse case when Ub is specified instead of U∞
# TODO: check if TU, TF should be used to enforce type compatibility in functions
# TODO: mention frame of reference for computeimpulse
# TODO: consider no deepcopy for new vortices in the methods and use deepcopy in the scripts instead
# TODO: redo PotentialFlow.jl solution in example 6 with uniform flow instead of moving plate

# TODO: create buffersol for computevortexvelocities?
# TODO: make computevortexvelocities in-place
# TODO: check computef̃limit for multiple bodies

"""
$(TYPEDEF)

Defines a grid-based vortex model with `Nb` bodies and `Ne` edges that are regularized. If `isshedding` is `true`, the strengths of the last `Ne` vortices in `vortices` can be computed in order to regularized the `Ne` edges.

# Fields

$(TYPEDFIELDS)

# Examples
(under construction)
"""
mutable struct VortexModel{Nb,Ne,TS<:AbstractPotentialFlowSystem}
    """g: The grid on which the vortex model is defined."""
    g::PhysicalGrid
    """bodies: Bodies in the vortex model."""
    bodies::StructVector{PotentialFlowBody}
    """vortices: Point vortices in the vortex model."""
    vortices::StructVector{Vortex}
    """U∞: Uniform flow in the vortex model."""
    U∞::Tuple{Float64,Float64}
    """system: Potential flow system that has to be solved with an `AbstractPotentialFlowRHS` and an `AbstractPotentialFlowSolution` to compute the potential flow that governs the vortex model.
    """
    system::TS

    """Internal fields"""
    _nodedata::Nodes{Dual}
    _edgedata::Edges{Primal}
    _bodyvectordata::VectorData
    _ψ::Nodes{Dual}
    _f::ScalarData
    _w::Nodes{Dual}
    _ψb::ScalarData
end

"""
$(TYPEDSIGNATURES)

Constructs a vortex model using the given function.
"""
function VortexModel(g::PhysicalGrid; bodies::StructVector{PotentialFlowBody}, vortices::StructVector{Vortex}=StructVector(Vortex[]), U∞::Tuple{Float64,Float64}=(0.0,0.0))

    vortices = deepcopy(vortices)

    e_idx = getregularizededges(bodies)

    Nb = length(bodies)
    Ne = length(e_idx)

    sizef = sum(length.(bodies))

    # Initialize data structures for internal use
    _nodedata = Nodes(Dual,size(g))
    _edgedata = Edges(Primal,size(g))
    _bodyvectordata = VectorData(sizef)
    _ψ = Nodes(Dual,size(g))
    _f = ScalarData(sizef)
    _w = Nodes(Dual,size(g))
    _ψb = ScalarData(sizef)

    L = plan_laplacian(size(_nodedata),with_inverse=true)
    Rbmat, Ebmat = _computeregularizationmatrix(g,VectorData(collect(bodies)), _f, _nodedata)

    one_vec = [ScalarData(sizef) for i in 1:Nb]
    for i in 1:Nb
        one_vec[i][getrange(bodies,i)] .= 1.0
    end

    if Ne == 0 # No regularized edges. Enforce circulation constraint.
        system = ConstrainedIBPoisson(L, Rbmat, Ebmat, one_vec, one_vec)
    else
        e_vec = [ScalarData(sizef) for i in 1:Ne]
        k = 0
        for i in 1:Nb
            for id in getregularizededges(bodies,i)
                k += 1
                e_vec[k][id] = 1.0
            end
        end
        if isempty(vortices) # With regularized edges, but no vortices. This is a steady case.
            system = SteadyRegularizedIBPoisson(L, Rbmat, Ebmat, one_vec, e_vec)
        else # With regularized edges and vortices. This is an unsteady case with vortex shedding.
            system = UnsteadyRegularizedIBPoisson(L, Rbmat, Ebmat, one_vec, e_vec)
        end
    end

    VortexModel{Nb,Ne,typeof(system)}(g, bodies, vortices, U∞, system, _nodedata, _edgedata, _bodyvectordata, _ψ, _f, _w, _ψb)
end

"""
$(SIGNATURES)

Replaces the vortices of the vortex model `vm` by a copy of `newvortices`.
"""
function setvortices!(vm::VortexModel{Nb,Ne}, newvortices::Vector{Vortex}) where {Nb,Ne}

    vm.vortices = StructArray(deepcopy(newvortices))
end

function setvortices!(vm::VortexModel{Nb,Ne}, newvortices::StructArray{Vortex}) where {Nb,Ne}

    vm.vortices = newvortices
end

"""
$(SIGNATURES)

Adds `newvortices` to the existing vortices in the vortex model `vm`.
"""
function pushvortices!(vm::VortexModel{Nb,Ne}, newvortices...) where {Nb,Ne}

    push!(vm.vortices, newvortices...)
end

"""
$(SIGNATURES)

Sets the positions of the vortices in the vortex model `vm` to the provided `X_vortices`.
"""
function setvortexpositions!(vm::VortexModel{Nb,Ne}, X_vortices::VectorData{Nv}) where {Nb,Ne,Nv}

    setpositions!(vm.vortices, X_vortices.u, X_vortices.v)
end

"""
$(SIGNATURES)

Returns the positions of the vortices in the vortex model `vm`.
"""
function getvortexpositions(vm::VortexModel{Nb,Ne}) where {Nb,Ne}

    return getpositions(vm.vortices)
end

function _updatesystemd_vec!(vm::VortexModel{Nb,Ne,UnsteadyRegularizedIBPoisson{Nb,Ne,TU,TF}}, idx::Vector{Int}) where {Nb,Ne,TU,TF}

    H = Regularize(getpositions(vm.vortices), cellsize(vm.g), I0=origin(vm.g), ddftype=CartesianGrids.M4prime, issymmetric=true)

    Γ = getstrengths(vm.vortices)
    for i in 1:Ne
        Γ .= 0
        Γ[idx[i]] = 1.0
        H(vm.system.d_vec[i],Γ)
    end
end

"""
$(SIGNATURES)

Returns the flow velocity, associated with the discrete vector potential field `ψ`, as `VectorData` at the locations of the vortices stored in the vortex model `vm`.
"""
function computevortexvelocities(vm::VortexModel{Nb,Ne}, ψ::Nodes{Dual}) where {Nb,Ne}

    H = Regularize(getpositions(vm.vortices), cellsize(vm.g), I0=origin(vm.g), ddftype=CartesianGrids.M4prime, issymmetric=true)

    Ẋ_vortices = VectorData(length(vm.vortices))
    # Velocity is the curl of the vector potential
    # The discrete curl operator requires dividing by the cellsize to account for the grid spacing
    curl!(_vm.edgedata, ψ)
    vm._edgedata ./= cellsize(vm.g)

    # For consistent interpolation, first interpolate the velocity to the nodes and use H to interpolate from the nodes to the vortices
    grid_interpolate!(_vm.nodedata, vm._edgedata.u);
    H(Ẋ_vortices.u, vm._nodedata)
    grid_interpolate!(vm._nodedata, vm._edgedata.v);
    H(Ẋ_vortices.v, vm._nodedata)

    return Ẋ_vortices
end

"""
$(SIGNATURES)

Returns the flow velocity as `VectorData` at the locations of the vortices stored in the vortex model `vm`, accounting for bodies in `vm` and conditions in `kwargs`. If the `vm` has `Ne` regularized edges and vortices, the strengths of the last `Ne` vortices will be computed and set in `vm`.
"""
function computevortexvelocities(vm::VortexModel{Nb,Ne}; kwargs...) where {Nb,Ne}

    # The strengths of the Ne last vortices will be calculated in solvesystem and should be set to zero before computing the vorticity field such that they are not included in w
    for k in 1:Ne
        vm.vortices.Γ[end-Ne+k] = 0.0
    end

    computew!(vm._w, vm)

    sol = solvesystem(vm, vm._w; kwargs...)

    Ẋ_vortices = computevortexvelocities(vm,sol.ψ)

    return Ẋ_vortices
end

"""
$(SIGNATURES)

Computes the vorticity field `w` associated with the vortices stored in the vortex model `vm` on the physical grid.
"""
function computew!(wphysical::Nodes{Dual}, vm::VortexModel{Nb,Ne})::Nodes{Dual} where {Nb,Ne}

    H = Regularize(getpositions(vm.vortices), cellsize(vm.g), I0=origin(vm.g), ddftype=CartesianGrids.M4prime, issymmetric=true)

    if isempty(vm.vortices)
        wphysical .= 0.0
        return wphysical
    end

    Γ = getstrengths(vm.vortices)
    H(wphysical,Γ)
    wphysical ./= cellsize(vm.g)^2 # Divide by the Δx² to ensure that ∫wdA = ΣΓ

    return wphysical
end

"""
$(SIGNATURES)

Computes the vorticity field associated with the vortices stored in the vortex model `vm` on the physical grid.
"""
function computew(vm::VortexModel{Nb,Ne})::Nodes{Dual} where {Nb,Ne}

    w = Nodes(Dual,size(vm.g))
    computew!(w,vm)

    return w
end

"""
$(SIGNATURES)

Computes the potential flow solution `sol` for the bodies and uniform flow in the vortex model `vm` and vorticity on the physical grid `wphysical`. If the vortex model contains bodies with regularized edges and a number of vortices which is greater than or equal to the number of regularized edges, the returned solution contains the computed strengths of the last N vortices in `vm`, with N the number of regularized edges.
"""
function solvesystem!(sol::ConstrainedIBPoissonSolution, vm::VortexModel{Nb,0,ConstrainedIBPoisson{Nb,TU,TF}}, wphysical::Nodes{Dual}) where {Nb,TU,TF}

    _computeψboundaryconditions!(vm._ψb, vm)
    vm._w .= wphysical
    Γb = deepcopy(vm.bodies.Γb)

    rhs = ConstrainedIBPoissonRHS(vm._w, vm._ψb, Γb)
    _scaletoindexspace!(rhs,cellsize(vm.g))
    ldiv!(sol,vm.system,rhs)
    _scaletophysicalspace!(sol,cellsize(vm.g))

    _addψ∞!(sol.ψ,vm) # Add the uniform flow to the approximation to the continuous stream function field

    return sol
end

function solvesystem!(sol::ConstrainedIBPoissonSolution, vm::VortexModel{Nb,Ne,SteadyRegularizedIBPoisson{Nb,Ne,TU,TF}}, wphysical::Nodes{Dual}) where {Nb,Ne,TU,TF}

    _computeψboundaryconditions!(vm._ψb, vm)

    vm._w .= wphysical

    f̃lim_vec = Vector{f̃Limit}()
    for i in 1:Nb
        append!(f̃lim_vec, _computef̃limit.(vm.bodies[i].σ, Ref(vm.bodies[i].points), sum(vm.system.f₀_vec[i])))
    end
    rhs = ConstrainedIBPoissonRHS(vm._w, vm._ψb, f̃lim_vec)
    _scaletoindexspace!(rhs, cellsize(vm.g))
    ldiv!(sol, vm.system, rhs)
    _scaletophysicalspace!(sol,cellsize(vm.g))

    _addψ∞!(sol.ψ,vm) # Add the uniform flow to the approximation to the continuous stream function field

    for i in 1:Nb
        vm.bodies[i].Γb = sum(view(sol.f,getrange(vm.bodies,i)))
    end

    return sol
end

function solvesystem!(sol::ConstrainedIBPoissonSolution, vm::VortexModel{Nb,Ne,UnsteadyRegularizedIBPoisson{Nb,Ne,TU,TF}}, wphysical::Nodes{Dual}) where {Nb,Ne,TU,TF}

    _computeψboundaryconditions!(vm._ψb, vm)
    _updatesystemd_vec!(vm,collect(length(vm.vortices)-(Ne-1):length(vm.vortices))) # Update the system.d_vec to agree with the latest vortex positions

    vm._w .= wphysical
    f̃lim_vec = Vector{f̃Limits}()
    for i in 1:Nb
        append!(f̃lim_vec, _computef̃limits.(vm.bodies[i].σ, Ref(vm.bodies[i].points), sum(vm.system.f₀_vec[i])))
    end
    Γw = -vm.bodies.Γb # Γw = ∫wdA

    rhs = UnsteadyRegularizedIBPoissonRHS(vm._w, vm._ψb, f̃lim_vec, Γw)
    _scaletoindexspace!(rhs,cellsize(vm.g))
    ldiv!(sol,vm.system,rhs)
    _scaletophysicalspace!(sol,cellsize(vm.g))

    _addψ∞!(sol.ψ,vm) # Add the uniform flow to the approximation to the continuous stream function field

    # set the strengths of the last Ne vortices
    for i in 1:Ne
        vm.vortices[end-Ne+i].Γ = sol.δΓ_vec[i]
    end
    subtractcirculation!(vm.bodies, sol.δΓ_vec)

    return sol
end

"""
$(SIGNATURES)

Computes and returns the potential flow solution for the bodies in the vortex model `vm` and vorticity on the physical grid `wphysical`. If the vm contains bodies with regularized edges and a number of vortices which is greater than or equal to the number of regularized edges, the returned solution contains the computed strengths of the last N vortices in `vm`, with N the number of regularized edges. A uniform flow, body velocities, bound circulation values for the unregularized bodies, and suction parameters for regularized bodies can be specified using the `parameters` keyword.
"""
function solvesystem(vm::VortexModel{Nb,0,ConstrainedIBPoisson{Nb,TU,TF}}, wphysical::Nodes{Dual}) where {Nb,TU,TF}

    sol = ConstrainedIBPoissonSolution(vm._ψ, vm._f, zeros(Float64,Nb), Float64[])
    solvesystem!(sol, vm, wphysical)

    return sol
end

function solvesystem(vm::VortexModel{Nb,Ne,SteadyRegularizedIBPoisson{Nb,Ne,TU,TF}}, wphysical::Nodes{Dual}) where {Nb,Ne,TU,TF}

    sol = ConstrainedIBPoissonSolution(vm._ψ, vm._f, zeros(Float64,Nb), Float64[])
    solvesystem!(sol, vm, wphysical)

    return sol
end

function solvesystem(vm::VortexModel{Nb,Ne,UnsteadyRegularizedIBPoisson{Nb,Ne,TU,TF}}, wphysical::Nodes{Dual}) where {Nb,Ne,TU,TF}

    sol = ConstrainedIBPoissonSolution(vm._ψ, vm._f, zeros(Float64,Nb), zeros(Float64,Nb))
    solvesystem!(sol, vm, wphysical)

    return sol
end

"""
$(SIGNATURES)

Computes and returns the stream function field on the physical grid for the potential flow associated with the current state of the vortex model `vm`.
"""
function computeψ(vm::VortexModel)

    computew!(vm._w, vm)
    sol = solvesystem(vm, vm._w)

    return sol.ψ
end

function _addψ∞!(ψ, vm::VortexModel)

    xg,yg = coordinates(ψ,vm.g)
    ψ .+= vm.U∞[1].*yg' .- vm.U∞[2].*xg

    return ψ
end

function _computeψboundaryconditions!(ψb::ScalarData, vm::VortexModel)

    computebodypointsvelocity!(vm._bodyvectordata.u, vm.bodies, 1) # Convert Ub into VectorData corresponding to the body points
    computebodypointsvelocity!(vm._bodyvectordata.v, vm.bodies, 2) # Convert Ub into VectorData corresponding to the body points

    # The discrete streamfunction field is constrained to a prescribed streamfunction on the body that describes the body motion. The body presence in the uniform flow is taken into account by subtracting its value from the body motion (i.e. a body motion in the -U∞ direction) and adding the uniform flow at the end of the solvesystem routine.
    ψb .= -vm.U∞[1]*(collect(vm.bodies)[2]) .+ vm.U∞[2]*(collect(vm.bodies)[1]);
    ψb .+= vm._bodyvectordata.u .* collect(vm.bodies)[2] .- vm._bodyvectordata.v .* collect(vm.bodies)[1]

    return ψb
end

function _computeψboundaryconditions!(ψb::ScalarData, vm::VortexModel{0})
    return
end

"""
$(SIGNATURES)

Computes the impulse associated with the vorticity `wphysical`, the bound vortex sheet strength `fphysical`, and the velocities `Ubvec` of the discrete points of the bodies in `vm`.
"""
function computeimpulse(vm::VortexModel{Nb,Ne}, wphysical::Nodes{Dual}, fphysical::ScalarData, Ubvec) where {Nb,Ne}

    @unpack g, vortices, bodies, _bodyvectordata = vm

    xg, yg = coordinates(wphysical,g)
    Δx = cellsize(g)

    # Formula 61 (see formula 6.16 in book)
    impulse = [Δx^2*sum(wphysical.*yg'),Δx^2*sum(-wphysical.*xg)]

    for i in 1:Nb
        impulse += _computeimpulsesurfaceintegral(bodies[i], fphysical[getrange(bodies,i)], Ubvec.u[getrange(bodies,i)], Ubvec.v[getrange(bodies,i)])
    end

    return impulse[1], impulse[2]
end

"""
$(SIGNATURES)

Computes the impulse associated with the current state vortices and bodies in the vortex model `vm`.
"""
function computeimpulse(vm::VortexModel)

    @unpack bodies, _nodedata, _bodyvectordata, _w = vm

    computew!(_nodedata,vm)
    # Solve system
    sol = solvesystem(vm, _nodedata, parameters=parameters)
    # Compute vorticity field again, because now it can contain the vorticity of newly shedded vortices
    computew!(_nodedata,vm)
    # Convert Ub into VectorData corresponding with the body points
    _computebodypointsvelocity!(_bodyvectordata,parameters.Ub,bodies)

    P_x, P_y = computeimpulse(vm, _nodedata, sol.f, _bodyvectordata)

    return P_x, P_y
end

function _computeimpulsesurfaceintegral(body::Body{N,RigidBodyTools.ClosedBody}, f, u, v) where {N}
    nx,ny = normalmid(body)
    Δs = dlengthmid(body)
    return _computecrossproductsurfaceintegral(body,(f./Δs + nx.*v - ny.*u).*Δs)
end

function _computeimpulsesurfaceintegral(body::Body{N,RigidBodyTools.OpenBody}, f, u, v) where {N}
    return _computecrossproductsurfaceintegral(body,f)
end

"""
$(SIGNATURES)

Computes the translational coefficients of the added mass matrix of the bodies in the vortex model `vm`.
"""
function computeaddedmassmatrix(vm::VortexModel{Nb,Ne}) where {Nb,Ne}
    @unpack g, vortices, bodies, _bodyvectordata, _w = vm

    # For now only translational
    M = zeros(Nb*2,Nb*2)

    computew!(_w,vm)

    for movingbodyindex in 1:Nb
        for dir in 1:2
            Ub = fill((0.0,0.0),Nb)
            if dir == 1
                Ub[movingbodyindex] = (1.0,0.0)
            else
                Ub[movingbodyindex] = (0.0,1.0)
            end
            _computebodypointsvelocity!(_bodyvectordata,Ub,bodies)
            sol = solvesystem(vm, _w)
            for i in 1:Nb
                M[(i-1)*2+1:(i-1)*2+2,(movingbodyindex-1)*2+dir] .= _computeimpulsesurfaceintegral(bodies[i], sol.f[getrange(bodies,i)], _bodyvectordata.u[getrange(bodies,i)], _bodyvectordata.v[getrange(bodies,i)])
            end
        end
    end
    return M
end

function show(io::IO, model::VortexModel{Nb,Ne,isshedding}) where {Nb,Ne,isshedding}
    NX = model.g.N[1]
    NY = model.g.N[2]
    N = length(model._f)
    Nv = length(model.vortices)

    println(io, "Vortex model on a grid of size $NX x $NY and $N immersed points with $((Nv == 1) ? "1 vortex" : "$Nv vortices"), $((Nb == 1) ? "1 body" : "$Nb bodies"), and $((Ne == 1) ? "1 regularized edge" : "$Ne regularized edges")")
end
