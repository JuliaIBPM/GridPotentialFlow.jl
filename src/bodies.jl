import Base: show

export PotentialFlowBody, subtractcirculation!, getΓ, setΓ!

"""
$TYPEDEF

Defines a `RigidBodyTools.Body` that combines a `RigidBodyTools.Body` with a circulation `Γ`, information about the edges `edges` and their corresponding suction parameters.
"""
mutable struct PotentialFlowBody{N,C,ST} <: Body{N,C} where {ST<:Union{Vector{SuctionParameter},Vector{SuctionParameterRange}}}
    """body: `RigidBodyTools.Body`"""
    body::Body{N,C}
    """Γ: Current circulation about the body."""
    Γ::Float64
    """regidx: Array of the indices of any regularized edges."""
    edges::Vector{Int}
    """σ: Array of suction parameters or suction parameter ranges associated with the entries of `edges`."""
    σ::ST
end

"""
$(TYPEDSIGNATURES)

Constructs a potential flow body from body `b` with initial circulation `Γ`. Edges that should be regularized can be specified in `edges` and their suction parameter ranges in `σ`, which defaults to zero for each regularized edge.
"""
function PotentialFlowBody(b::Body{N,C}; Γ = 0.0, edges = Int64[], σ = SuctionParameterRange[]) where {N,C}

    Ne = length(edges)

    if Ne > 0 && isempty(σ)
        σ = [SuctionParameterRange(0.0,0.0) for e in 1:Ne]
    elseif Ne > length(σ)
        println("Not enough suction parameters provided. Setting all suction parameters to zero.")
    end

    return PotentialFlowBody(b,Γ,edges,σ)
end

for f in [:diff,:endpoints,:midpoints,:centraldiff]
    _f = Symbol("_"*string(f))
    @eval RigidBodyTools.$f(pfb::PotentialFlowBody{N,C,ST};ref=false) where {N,C<:RigidBodyTools.BodyClosureType,ST} = RigidBodyTools.$_f(pfb.body,Val(ref))
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

Sets the bound circulation of the body `b` to `Γ`.
"""
function setΓ!(b::PotentialFlowBody, Γ::Float64)
    b.Γ = Γ
end

# """
# $(TYPEDSIGNATURES)
#
# Returns the indices in the global set of surface point data of the regularized points of body `i` in `b`.
# """
# function getregularizededges(b::AbstractVector{PotentialFlowBody},i::Int)
#     r = getrange(b,i)
#     idx = deepcopy(b[i].edges) .+ r[1] .- 1
#     return idx
# end
#
# """
# $(TYPEDSIGNATURES)
#
# Returns the indices in the global set of surface point data of all regularized points in `b`.
# """
# function getregularizededges(b::AbstractVector{PotentialFlowBody})
#     idx = Int64[]
#     for i in 1:length(b)
#         append!(idx,getregularizededges(b,i))
#     end
#     return idx
# end
#
# """
# $(TYPEDSIGNATURES)
#
#
# """
# function subtractcirculation!(b::AbstractVector{PotentialFlowBody}, δΓ_vec::AbstractVector)
#     i = 0
#     for j in 1:length(b)
#         for k in 1:length(b[j].edges)
#             i += 1
#             b[j].Γ -= δΓ_vec[i]
#         end
#     end
# end



function Base.show(io::IO, b::PotentialFlowBody)
    Ne = length(b.edges)
    iobody = IOBuffer()
    Base.show(iobody, b.body)
    sbody = String(take!(iobody))
    println(io,"Potential flow body with $Ne regularized $(Ne == 1 ? "edge" : "edges")")
    println(io,"Shape and position: $(sbody[1:end-1])")
    println(io,"Current circulation: $(b.Γ)")
    if Ne > 0
        println(io,"Suction parameter limits:")
        for e in 1:Ne
            println(io,"   Edge $e: $(b.σ[e])")
        end
    end
end

@recipe f(::Type{PotentialFlowBody}, pfb::PotentialFlowBody) = pfb.points
