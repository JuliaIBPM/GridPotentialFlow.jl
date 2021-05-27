import CartesianGrids: curl!, Laplacian
import LinearAlgebra: Diagonal, norm
import Base: show

using RecipesBase

export VortexModel, streamfunction, vorticity, vorticity!, vortexvelocities!, _computeregularizationmatrix, setvortexstrengths!, setvortexpositions!, getvortexpositions, setvortices!, pushvortices!, impulse, addedmass, solve, solve!

# TODO: check if _streamfunctionbcs needs to be faster
# TODO: check impulse cases for both Ub and U∞ specified
# TODO: check if TU, TF should be used to enforce type compatibility in functions
# TODO: mention frame of reference for computeimpulse
# TODO: consider no deepcopy for new vortices in the methods and use deepcopy in the scripts instead
# TODO: redo PotentialFlow.jl solution in example 6 with uniform flow instead of moving plate

# TODO: check computef̃limit for multiple bodies (f₀). Maybe Γ₀ needs to be a vector.
# TODO: check if systems can be reduced to basic one and other one for all the relevant cases.
# TODO: set vortex strengths in script? Then warn in impulse that it does not regularize
# TODO: case where steady regularized bodies are next to unregularized bodies with circulation.
# TODO: change vortexvelocities to accept w and don't set vortices to 0 anymore.
# TODO: add parameter to VortexModel constructor to make model unsteady even if no vortices are provided.

"""
$(TYPEDEF)

Defines a grid-based vortex model with `Nb` bodies and `Ne` edges that are regularized. If `isshedding` is `true`, the strengths of the last `Ne` vortices in `vortices` can be computed in order to regularized the `Ne` edges.

# Fields

$(TYPEDFIELDS)

# Examples
(under construction)
"""
mutable struct VortexModel{Nb,Ne,TS<:Union{AbstractPotentialFlowSystem,Laplacian},TU<:Nodes,TE<:Edges,TF<:ScalarData,TX<:VectorData}
    """g: The grid on which the vortex model is defined."""
    g::PhysicalGrid
    """bodies: Bodies in the vortex model."""
    bodies::Vector{PotentialFlowBody}
    """vortices: Point vortices in the vortex model."""
    vortices::StructVector{Vortex}
    """U∞: Uniform flow in the vortex model."""
    U∞::Tuple{Float64,Float64}
    """system: Potential flow system that has to be solved with an `AbstractPotentialFlowRHS` and an `AbstractPotentialFlowSolution` to compute the potential flow that governs the vortex model.
    """
    system::TS

    """Internal fields"""
    _nodedata::TU
    _edgedata::TE
    _bodyvectordata::TX
    _ψ::TU
    _f::TF
    _w::TU
    _ψb::TF
end

"""
$(TYPEDSIGNATURES)

Constructs a vortex model using the given function.
"""
function VortexModel(g::PhysicalGrid; bodies::Vector{PotentialFlowBody}=Vector{PotentialFlowBody}(), vortices::StructVector{Vortex}=StructVector(Vortex[]), U∞::Tuple{Float64,Float64}=(0.0,0.0))

    vortices = deepcopy(vortices)

    e_idx = getregularizededges(bodies)

    Nb = length(bodies)
    Ne = length(e_idx)

    if isempty(bodies)
        sizef = 0
    else
        sizef = sum(length.(bodies))
    end

    # Initialize data structures for internal use
    _nodedata = Nodes(Dual,size(g))
    _edgedata = Edges(Primal,size(g))
    _bodyvectordata = VectorData(sizef)
    _ψ = Nodes(Dual,size(g))
    _f = ScalarData(sizef)
    _w = Nodes(Dual,size(g))
    _ψb = ScalarData(sizef)

    L = plan_laplacian(size(_nodedata),with_inverse=true)
    if Nb == 0
        system = L
    else
        regop = Regularize(VectorData(collect(bodies)), cellsize(g), I0=origin(g), ddftype = CartesianGrids.Yang3, issymmetric=true)
        Rmat,_ = RegularizationMatrix(regop, _f, _nodedata)
        Emat = InterpolationMatrix(regop, _nodedata, _f)

        one_vec = [ScalarData(sizef) for i in 1:Nb]
        for i in 1:Nb
            one_vec[i][getrange(bodies,i)] .= 1.0
        end

        if Ne == 0 # System without regularized edges. Enforce circulation constraints.
            system = ConstrainedIBPoisson(L, Rmat, Emat, one_vec, one_vec)
        else # System with regularized edges. Enforce edge constraints.
            e_vec = [ScalarData(sizef) for i in 1:Ne]
            k = 0
            for i in 1:Nb
                for id in getregularizededges(bodies,i)
                    k += 1
                    e_vec[k][id] = 1.0
                end
            end
            if isempty(vortices) # With regularized edges, but no vortices. This is a steady case.
                system = SteadyRegularizedIBPoisson(L, Rmat, Emat, one_vec, e_vec)
            else # With regularized edges and vortices. This is an unsteady case with vortex shedding.
                vidx_vec = Vector{Vector{Int64}}()
                vidx = collect(1:Ne)
                for i in 1:Nb
                    vidx_i = [pop!(vidx) for k in 1:length(bodies[i].edges)]
                    push!(vidx_vec, vidx_i)
                end
                system = UnsteadyRegularizedIBPoisson(L, Rmat, Emat, one_vec, e_vec, vidx_vec)
            end
        end
    end

    VortexModel{Nb,Ne,typeof(system),typeof(_ψ),typeof(_edgedata),typeof(_f),typeof(_bodyvectordata)}(g, bodies, vortices, U∞, system, _nodedata, _edgedata, _bodyvectordata, _ψ, _f, _w, _ψb)
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

Sets the positions of the vortices in the vortex model `vm` to the provided `X`.
"""
function setvortexpositions!(vm::VortexModel{Nb,Ne}, X::VectorData{Nv}) where {Nb,Ne,Nv}

    setpositions!(vm.vortices, X.u, X.v)
end

"""
$(SIGNATURES)

Sets the strenghts of the vortices in the vortex model `vm` to the provided `Γ`. Indices can be provided in `idx` to specify the indices of the vortices whose strengths have to be set, in which case the length of `Γ` must be the same as the number of indices.
"""
function setvortexstrengths!(vm::VortexModel{Nb,Ne}, Γ::TΓ, idx=nothing) where {Nb,Ne,TΓ<:AbstractVector}

    setstrengths!(vm.vortices, Γ, idx)
end

"""
$(SIGNATURES)

Returns the positions of the vortices in the vortex model `vm`.
"""
function getvortexpositions(vm::VortexModel{Nb,Ne}) where {Nb,Ne}

    return getpositions(vm.vortices)
end

function _updatesystemd_vec!(vm::VortexModel{Nb,Ne,UnsteadyRegularizedIBPoisson{Nb,Ne,TU,TF}}, idx::Vector{Int}) where {Nb,Ne,TU,TF}

    v = view(vm.vortices,idx)

    H = Regularize(getpositions(v), cellsize(vm.g), I0=origin(vm.g), ddftype=CartesianGrids.M4prime, issymmetric=true)

    Γ = getstrengths(v)
    for i in 1:Ne
        Γ .= 0
        Γ[i] = 1.0
        H(vm.system.d_vec[i],Γ)
    end
end

"""
$(SIGNATURES)

Computes and stores in `Ẋ_vortices` the flow velocity, associated with the discrete vector potential field `ψ`, as `VectorData` at the locations of the vortices stored in the vortex model `vm`.
"""
function vortexvelocities!(Ẋ_vortices::VectorData, vm::VortexModel{Nb,Ne,TS,TU}, ψ::TU) where {Nb,Ne,TS,TU}

    H = Regularize(getpositions(vm.vortices), cellsize(vm.g), I0=origin(vm.g), ddftype=CartesianGrids.M4prime, issymmetric=true)

    # Velocity is the curl of the vector potential
    # The discrete curl operator requires dividing by the cellsize to account for the grid spacing
    curl!(vm._edgedata, ψ)
    vm._edgedata ./= cellsize(vm.g)

    # For consistent interpolation, first interpolate the velocity to the nodes and use H to interpolate from the nodes to the vortices
    grid_interpolate!(vm._nodedata, vm._edgedata.u);
    H(Ẋ_vortices.u, vm._nodedata)
    grid_interpolate!(vm._nodedata, vm._edgedata.v);
    H(Ẋ_vortices.v, vm._nodedata)

    return Ẋ_vortices
end

"""
$(SIGNATURES)

Returns the flow velocity as `VectorData` at the locations of the vortices stored in the vortex model `vm`, accounting for bodies in `vm`. If the `vm` has `Ne` regularized edges and vortices, the strengths of the last `Ne` vortices will be computed and set in `vm`.
"""
function vortexvelocities!(vm::VortexModel{Nb,Ne}) where {Nb,Ne}

    # The strengths of the Ne last vortices will be calculated in solve and should be set to zero before computing the vorticity field such that they are not included in w
    for k in 1:Ne
        vm.vortices.Γ[end-Ne+k] = 0.0
    end

    vorticity!(vm._w, vm)
    sol = solve(vm, vm._w)
    Ẋ_vortices = VectorData(length(vm.vortices))

    for k in 1:Ne
        vm.vortices.Γ[end-Ne+k] = sol.δΓ_vec[k]
    end

    return vortexvelocities!(Ẋ_vortices, vm, sol.ψ)
end

"""
$(SIGNATURES)

Computes the vorticity field `w` associated with the vortices stored in the vortex model `vm` on the physical grid.
"""
function vorticity!(wphysical::TU, vm::VortexModel{Nb,Ne,TS,TU}) where {Nb,Ne,TS,TU}

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
function vorticity(vm::VortexModel{Nb,Ne}) where {Nb,Ne}

    w = Nodes(Dual,size(vm.g))
    vorticity!(w,vm)

    return w
end

"""
$(SIGNATURES)

Computes the potential flow solution `sol` for the bodies and uniform flow in the vortex model `vm` and vorticity on the physical grid `wphysical`. If the vortex model contains bodies with regularized edges and a number of vortices which is greater than or equal to the number of regularized edges, the returned solution contains the computed strengths of the last N vortices in `vm`, with N the number of regularized edges.
"""
function solve!(sol::PoissonSolution, vm::VortexModel{0,0,TS,TU}, wphysical::TU) where {Nb,TU,TS<:Laplacian}

    vm._w .= wphysical
    rhs = PoissonRHS(vm._w)
    _scaletoindexspace!(rhs, cellsize(vm.g))
    ldiv!(sol.ψ, vm.system, rhs.w)
    _scaletophysicalspace!(sol, cellsize(vm.g))

    _addψ∞!(sol.ψ, vm) # Add the uniform flow to the approximation to the continuous stream function field

    return sol
end

function solve!(sol::ConstrainedIBPoissonSolution, vm::VortexModel{Nb,0,ConstrainedIBPoisson{Nb,TU,TF}}, wphysical::TU) where {Nb,TU,TF}

    _streamfunctionbcs!(vm._ψb, vm)
    vm._w .= wphysical
    Γb = deepcopy(getΓb(vm.bodies))

    rhs = ConstrainedIBPoissonRHS(vm._w, vm._ψb, Γb)
    _scaletoindexspace!(rhs, cellsize(vm.g))
    ldiv!(sol, vm.system, rhs)
    _scaletophysicalspace!(sol, cellsize(vm.g))

    _addψ∞!(sol.ψ, vm) # Add the uniform flow to the approximation to the continuous stream function field

    return sol
end

function solve!(sol::ConstrainedIBPoissonSolution, vm::VortexModel{Nb,Ne,SteadyRegularizedIBPoisson{Nb,Ne,TU,TF}}, wphysical::TU) where {Nb,Ne,TU,TF}

    _streamfunctionbcs!(vm._ψb, vm)

    vm._w .= wphysical

    f̃lim_vec = Vector{f̃Limit}()
    for i in 1:Nb
        append!(f̃lim_vec, _computef̃limit.(vm.bodies[i].σ, Ref(vm.bodies[i].points), sum(vm.system.f₀_vec[i])))
    end
    rhs = ConstrainedIBPoissonRHS(vm._w, vm._ψb, f̃lim_vec)
    _scaletoindexspace!(rhs, cellsize(vm.g))
    ldiv!(sol, vm.system, rhs)
    _scaletophysicalspace!(sol, cellsize(vm.g))

    _addψ∞!(sol.ψ, vm) # Add the uniform flow to the approximation to the continuous stream function field

    # This next for loop should be done in the script
    # for i in 1:Nb
    #     vm.bodies[i].Γb = sum(view(sol.f,getrange(vm.bodies,i)))
    # end

    return sol
end

function solve!(sol::ConstrainedIBPoissonSolution, vm::VortexModel{Nb,Ne,UnsteadyRegularizedIBPoisson{Nb,Ne,TU,TF}}, wphysical::TU) where {Nb,Ne,TU,TF}

    _streamfunctionbcs!(vm._ψb, vm)
    _updatesystemd_vec!(vm,collect(length(vm.vortices)-(Ne-1):length(vm.vortices))) # Update the system.d_vec to agree with the latest vortex positions

    vm._w .= wphysical
    f̃lim_vec = Vector{f̃Limits}()
    for i in 1:Nb
        append!(f̃lim_vec, _computef̃limits.(vm.bodies[i].σ, Ref(vm.bodies[i].points), sum(vm.system.f₀_vec[i])))
    end
    Γw = -getΓb(vm.bodies) # Γw = ∫wdA

    rhs = UnsteadyRegularizedIBPoissonRHS(vm._w, vm._ψb, f̃lim_vec, Γw)
    _scaletoindexspace!(rhs, cellsize(vm.g))
    ldiv!(sol, vm.system, rhs)
    _scaletophysicalspace!(sol, cellsize(vm.g))

    _addψ∞!(sol.ψ, vm) # Add the uniform flow to the approximation to the continuous stream function field

    return sol
end

"""
$(SIGNATURES)

Computes and returns the potential flow solution for the bodies in the vortex model `vm` and vorticity on the physical grid `wphysical`. If the vm contains bodies with regularized edges and a number of vortices which is greater than or equal to the number of regularized edges, the returned solution contains the computed strengths of the last N vortices in `vm`, with N the number of regularized edges.
"""
function solve(vm::VortexModel{0,0,TS,TU}, wphysical::TU) where {Nb,Ne,TS<:Laplacian,TU}

    sol = PoissonSolution(TU())
    solve!(sol, vm, wphysical)

    return sol
end

function solve(vm::VortexModel{Nb,Ne,TS,TU,TE,TF}, wphysical::TU) where {Nb,Ne,TS<:AbstractPotentialFlowSystem,TU,TE,TF}

    sol = ConstrainedIBPoissonSolution(TU(), TF(), zeros(Float64,Nb), zeros(Float64,Ne))
    solve!(sol, vm, wphysical)

    return sol
end

"""
$(SIGNATURES)

Computes the stream function field `ψ` on the physical grid for the potential flow associated with the current state of the vortex model `vm` and vorticity on the physical grid `wphysical`.
"""
function streamfunction!(ψ::TU, vm::VortexModel{0,0,TS,TU}, wphysical::TU) where {Nb,Ne,TS<:Laplacian,TU}

    sol = PoissonSolution(ψ)
    solve!(sol, vm, wphysical)

    return sol.ψ
end

function streamfunction!(ψ::TU, vm::VortexModel{Nb,Ne,TS,TU,TE,TF}, wphysical::TU) where {Nb,Ne,TS<:AbstractPotentialFlowSystem,TU,TE,TF}

    sol = ConstrainedIBPoissonSolution(ψ, _f, zeros(Float64,Nb), zeros(Float64,Ne))
    solve!(sol, vm, wphysical)

    return sol.ψ
end


"""
$(SIGNATURES)

Computes and returns the stream function field on the physical grid for the potential flow associated with the current state of the vortex model `vm`.
"""
function streamfunction(vm::VortexModel)

    vorticity!(vm._w, vm)
    sol = solve(vm, vm._w)

    return sol.ψ
end

"""
$(SIGNATURES)

Computes the impulse associated with the vorticity `wphysical`, the bound vortex sheet strength `fphysical`, and the velocities `Ubvec` of the discrete points of the bodies in `vm`.
"""
function impulse(vm::VortexModel{Nb,Ne}, wphysical::Nodes{Dual}, fphysical::ScalarData) where {Nb,Ne}

    xg, yg = coordinates(wphysical, vm.g)
    Δx = cellsize(vm.g)
    _bodypointsvelocity!(vm._bodyvectordata, vm.bodies)

    # Formula 61 (see formula 6.16 in book)
    p = [Δx^2*sum(wphysical.*yg'),Δx^2*sum(-wphysical.*xg)]

    for i in 1:Nb
        r = getrange(vm.bodies,i)
        p += _impulsesurfaceintegral(vm.bodies[i].points, fphysical[r], vm._bodyvectordata.u[r], vm._bodyvectordata.v[r])
    end

    return p[1], p[2]
end

"""
$(SIGNATURES)

Computes the impulse associated with the body motions in the vortex model `vm` and the vortex sheet strength in `sol`.
"""
function impulse(sol::TS, vm::VortexModel) where {TS<:AbstractIBPoissonSolution}
    # Compute vorticity field again, because now it can contain the vorticity of newly shedded vortices
    vorticity!(vm._nodedata, vm)
    impulse(vm, vm._nodedata, sol.f)
end

"""
$(SIGNATURES)

Computes the impulse associated with the current state vortices and bodies in the vortex model `vm`.
"""
function impulse(vm::VortexModel)
    vorticity!(vm._nodedata, vm)
    sol = solve(vm, vm._nodedata)
    impulse(sol, vm)
end

"""
$(SIGNATURES)

Computes the translational coefficients of the added mass matrix of the bodies in the vortex model `vm`.
"""
function addedmass(vm::VortexModel{Nb,Ne}) where {Nb,Ne}

    # Deepcopy bodies because we will change the velocities one by one
    oldbodies = deepcopy(vm.bodies)

    # For now only translational
    M = zeros(Nb*2,Nb*2)

    vorticity!(vm._w, vm)

    for movingbodyindex in 1:Nb
        for dir in 1:2
            for i in Nb
                vm.bodies[i].Ub = (0.0,0.0)
            end
            if dir == 1
                vm.bodies[movingbodyindex].Ub = (1.0,0.0)
            else
                vm.bodies[movingbodyindex].Ub = (0.0,1.0)
            end
            _bodypointsvelocity!(vm._bodyvectordata, vm.bodies)
            sol = solve(vm, vm._w)
            for i in 1:Nb
                r = getrange(vm.bodies,i)
                M[(i-1)*2+1:(i-1)*2+2,(movingbodyindex-1)*2+dir] .= _impulsesurfaceintegral(vm.bodies[i].points, sol.f[r], vm._bodyvectordata.u[r], vm._bodyvectordata.v[r])
            end
        end
    end
    vm.bodies = oldbodies
    return M
end

function _impulsesurfaceintegral(body::Body{N,RigidBodyTools.ClosedBody}, f, u, v) where {N}
    nx,ny = normalmid(body)
    Δs = dlengthmid(body)
    return _crossproductsurfaceintegral(body,(f./Δs + nx.*v - ny.*u).*Δs)
end

function _impulsesurfaceintegral(body::Body{N,RigidBodyTools.OpenBody}, f, u, v) where {N}
    return _crossproductsurfaceintegral(body,f)
end

function _crossproductsurfaceintegral(body::Body, z)
    rx,ry = collect(body)
    return [ry'*z, -rx'*z]
end

function _bodypointsvelocity!(v::VectorData, bodies)
    bodypointsvelocity!(v.u, bodies, 1)
    bodypointsvelocity!(v.v, bodies, 2)
    return v
end


function _addψ∞!(ψ, vm::VortexModel)

    xg,yg = coordinates(ψ,vm.g)
    ψ .+= vm.U∞[1].*yg' .- vm.U∞[2].*xg

    return ψ
end

function _streamfunctionbcs!(ψb::ScalarData, vm::VortexModel)

    _bodypointsvelocity!(vm._bodyvectordata, vm.bodies) # Convert Ub into VectorData corresponding to the body points

    # The discrete streamfunction field is constrained to a prescribed streamfunction on the body that describes the body motion. The body presence in the uniform flow is taken into account by subtracting its value from the body motion (i.e. a body motion in the -U∞ direction) and adding the uniform flow at the end of the solve routine.
    ψb .= -vm.U∞[1]*(collect(vm.bodies)[2]) .+ vm.U∞[2]*(collect(vm.bodies)[1]);
    ψb .+= vm._bodyvectordata.u .* collect(vm.bodies)[2] .- vm._bodyvectordata.v .* collect(vm.bodies)[1]

    return ψb
end

function _streamfunctionbcs!(ψb::ScalarData, vm::VortexModel{0})
    return
end

function show(io::IO, model::VortexModel{Nb,Ne,isshedding}) where {Nb,Ne,isshedding}
    NX = model.g.N[1]
    NY = model.g.N[2]
    N = length(model._f)
    Nv = length(model.vortices)

    println(io, "Vortex model on a grid of size $NX x $NY and $N immersed points with $((Nv == 1) ? "1 vortex" : "$Nv vortices"), $((Nb == 1) ? "1 body" : "$Nb bodies"), and $((Ne == 1) ? "1 regularized edge" : "$Ne regularized edges")")
end

@recipe function f(h::VortexModel)
    xlims := h.g.xlim[1]
    ylims := h.g.xlim[2]
    legend := :false
    grid := :false
    aspect_ratio := 1

    @series begin
        seriestype := :scatter
        markersize := 3
        markerstrokewidth := 0
        marker_z := h.vortices.Γ
        h.vortices.x, h.vortices.y
    end

    for b in h.bodies
        @series begin
            b.points
        end
    end
end
