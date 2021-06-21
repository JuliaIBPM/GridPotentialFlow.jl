import Base: show, length, collect
import RigidBodyTools: RigidTransform

export PotentialFlowBody, subtractcirculation!, getU, getΩ, getΓ, setU, setΩ, setΓ

"""
$TYPEDEF

Defines a rigid body for `VortexModel` in `GridPotentialFlow.jl`.
"""
mutable struct PotentialFlowBody
    """body: """
    points::Body
    """U: Linear velocity of the body."""
    U::Tuple{Float64,Float64}
    """Ω: Rotational velocity of the body."""
    Ω::Float64
    """Γ: Current circulation about the body."""
    Γ::Float64
    """regidx: Array of the indices of any regularized edges."""
    edges::Vector{Int}
    """σ: Array of suction parameters or suction parameter ranges associated with the entries of `edges`."""
    σ::Vector{Union{SuctionParameter,SuctionParameterRange}}
end

"""
$(TYPEDSIGNATURES)

Constructs a potential flow body with shape `b`, linear velocity `U`, rotational velocity `Ω`, and initial circulation `Γ`. Edges that should be regularized can be specified in `edges` and their suction parameter ranges in `σ`, which defaults to zero for each regularized edge.
"""
function PotentialFlowBody(b::Body{N,C}; U = (0.0,0.0), Ω = 0.0, Γ = 0.0, edges = Int64[], σ = SuctionParameterRange[]) where {N,C}

    Ne = length(edges)

    if Ne > 0 && isempty(σ)
        σ = [SuctionParameterRange(0.0,0.0) for e in 1:Ne]
    elseif Ne > length(σ)
        println("Not enough suction parameters provided. Setting all suction parameters to zero.")
    end

    return PotentialFlowBody(b,U,Ω,Γ,edges,σ)
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

Returns the linear velocity of the body `b`.
"""
function getU(b::PotentialFlowBody)
    return b.U
end

"""
$(TYPEDSIGNATURES)

Returns the rotational velocity of the body `b`.
"""
function getΩ(b::PotentialFlowBody)
    return b.Ω
end

"""
$(TYPEDSIGNATURES)

Returns the bound circulation of the body `b`.
"""
function getΓ(b::PotentialFlowBody)
    return b.Γ
end

"""
$(TYPEDSIGNATURES)

Sets the linear velocity of the body `b` to `U`.
"""
function setU(b::PotentialFlowBody, U::Tuple{Float64,Float64})
    b.U = U
end

"""
$(TYPEDSIGNATURES)

Sets the rotational velocity of the body `b` to `Ω`.
"""
function setΩ(b::PotentialFlowBody, Ω::Float64)
    b.Ω = Ω
end

"""
$(TYPEDSIGNATURES)

Sets the bound circulation of the body `b` to `Γ`.
"""
function setΓ(b::PotentialFlowBody, Γ::Float64)
    b.Γ = Γ
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
            b[j].Γ -= δΓ_vec[i]
        end
    end
end

function bodypointsvelocity!(v::AbstractVector, b::AbstractVector{PotentialFlowBody}, dir::Int)
    for i in 1:length(b)
        vi = view(v,getrange(b,i))
        vi .= b[i].U[dir]
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

function (T::RigidBodyTools.RigidTransform)(b::PotentialFlowBody)
    T(b.points)
    return b
end

@inline RigidBodyTools.dlength(b::PotentialFlowBody) = dlength(b.points)

@inline RigidBodyTools.dlengthmid(b::PotentialFlowBody) = dlengthmid(b.points)

@inline RigidBodyTools.normalmid(b::PotentialFlowBody) = normalmid(b.points)

function RigidBodyTools.dlengthmid(bl::AbstractVector{T}) where T<:PotentialFlowBody
    ds = Float64[]
    for b in bl
        dsb = dlengthmid(b)
        append!(ds,dsb)
    end
    return ds
end

function RigidBodyTools.normalmid(bl::AbstractVector{T}) where T<:PotentialFlowBody
    nx = Float64[]
    ny = Float64[]
    for b in bl
        nxb, nyb = normalmid(b)
        append!(nx,nxb)
        append!(ny,nyb)
    end
    return nx, ny
end

ImmersedLayers.areas(b::PotentialFlowBody) = areas(b.points)
ImmersedLayers.normals(b::PotentialFlowBody) = normals(b.points)

ImmersedLayers.areas(bl::AbstractVector{T}) where T<:PotentialFlowBody = dlengthmid(bl)
ImmersedLayers.normals(bl::AbstractVector{T}) where T<:PotentialFlowBody = VectorData(normalmid(bl))


function Base.show(io::IO, b::PotentialFlowBody)
    Ne = length(b.edges)
    iobody = IOBuffer()
    Base.show(iobody, b.points)
    sbody = String(take!(iobody))
    println(io,"Potential flow body with $Ne regularized $(Ne == 1 ? "edge" : "edges")")
    println(io,"Shape and position: $(sbody[1:end-1])")
    println(io,"Linear velocity: $(b.U)")
    println(io,"Angular velocity: $(b.Ω)")
    println(io,"Current circulation: $(b.Γ)")
    if Ne > 0
        println(io,"Suction parameter limits:")
        for e in 1:Ne
            println(io,"   Edge $e: $(b.σ[e])")
        end
    end
end

@recipe f(::Type{PotentialFlowBody}, pfb::PotentialFlowBody) = pfb.points
