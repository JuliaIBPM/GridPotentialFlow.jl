import Base: show, length, collect

export PotentialFlowBody

"""
$TYPEDEF

Defines a rigid body for `VortexModel` in `GridPotentialFlow.jl`.
"""
mutable struct PotentialFlowBody
    """body: """
    points::Body
    """Ub: Linear velocity of the body."""
    Ub::Tuple{Float64,Float64}
    """Ωb: Rotational velocity of the body."""
    Ωb::Float64
    """Γb: Current circulation about the body."""
    Γb::Float64
    """regidx: Array of the indices of any regularized edges."""
    edges::Vector{Int}
    """σ: Array of suction parameters or suction parameter ranges associated with the entries of `edges`."""
    σ::Vector{Union{SuctionParameter,SuctionParameterRange}}
end

"""
$(TYPEDSIGNATURES)

Constructs a potential flow body with shape `b`, linear velocity `Ub`, rotational velocity `Ωb`, and initial circulation `Γb`. Edges that should be regularized can be specified in `edges` and their suction parameter ranges in `σ`, which defaults to zero for each regularized edge.
"""
function PotentialFlowBody(b::Body{N,C}; Ub = (0.0,0.0), Ωb = 0.0, Γb = 0.0, edges = Int64[], σ = SuctionParameterRange[]) where {N,C}

    Ne = length(edges)

    if Ne > 0 && isempty(σ)
        σ = [SuctionParameterRange(0.0,0.0) for e in 1:Ne]
    elseif Ne > length(σ)
        println("Not enough suction parameters provided. Setting all suction parameters to zero.")
    end

    return PotentialFlowBody(b,Ub,Ωb,Γb,edges,σ)
end

"""
$(TYPEDSIGNATURES)

Returns the number of surface points in the body `b`.
"""
function Base.length(b::PotentialFlowBody)
    return length(b.points)
end

"""
$(TYPEDSIGNATURES)

Returns the subrange of indices in the global set of surface point data corresponding to body `i` in `b`.
"""
function getrange(b::AbstractVector{PotentialFlowBody},i::Int)
    first = 1
    j = 1
    while j < i
        first += length(b[j])
        j += 1
    end
    last = first+length(b[i])-1
    return first:last
end

"""
$(TYPEDSIGNATURES)

Returns a vector of the bound circulation values of the bodies in `b`.
"""
function getΓb(b::AbstractVector{PotentialFlowBody})
    return map(x->x.Γb,b)
end

"""
$(TYPEDSIGNATURES)

Returns the indices in the global set of surface point data of the regularized points of body `i` in `b`.
"""
function getregularizededges(b::AbstractVector{PotentialFlowBody},i::Int)
    r = getrange(b,i)
    idx = deepcopy(b[i].edges) .+ r[1] .- 1
    return idx
end

"""
$(TYPEDSIGNATURES)

Returns the indices in the global set of surface point data of all regularized points in `b`.
"""
function getregularizededges(b::AbstractVector{PotentialFlowBody})
    idx = Int64[]
    for i in 1:length(b)
        append!(idx,getregularizededges(b,i))
    end
    return idx
end

"""
$(TYPEDSIGNATURES)

Returns the indices in the global set of surface point data of all regularized points in `b`.
"""
function subtractcirculation!(b::AbstractVector{PotentialFlowBody}, δΓ_vec::AbstractVector)
    i = 0
    for j in 1:length(b)
        for k in 1:length(b[j].edges)
            i += 1
            b[j].Γb -= δΓ_vec[i]
        end
    end
end

function computebodypointsvelocity!(v::AbstractVector, b::AbstractVector{PotentialFlowBody}, dir::Int)
    for i in 1:length(b)
        vi = view(v,getrange(b,i))
        vi .= b[i].Ub[dir]
    end
end

function Base.collect(b::AbstractVector{PotentialFlowBody})
    x = Vector{Float64}()
    y = Vector{Float64}()
    for i in 1:length(b)
        xi, yi = RigidBodyTools.collect(b[i].points)
        append!(x,deepcopy(xi))
        append!(y,deepcopy(yi))
    end

    return x, y
end

function Base.show(io::IO, b::PotentialFlowBody)
    Ne = length(b.edges)
    iobody = IOBuffer()
    Base.show(iobody, b.points)
    sbody = String(take!(iobody))
    println(io,"Potential flow body with $Ne regularized $(Ne == 1 ? "edge" : "edges")")
    println(io,"Shape and position: $(sbody[1:end-1])")
    println(io,"Linear velocity: $(b.Ub)")
    println(io,"Angular velocity: $(b.Ωb)")
    println(io,"Current circulation: $(b.Γb)")
    if Ne > 0
        println(io,"Suction parameter limits:")
        for e in 1:Ne
            println(io,"   Edge $e: $(b.σ[e])")
        end
    end
end
