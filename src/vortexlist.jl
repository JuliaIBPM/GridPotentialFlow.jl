using RecipesBase

export VortexList

import Base: @propagate_inbounds,getindex, setindex!,iterate,size,length,push!,
              collect, view, vcat, lastindex

"""
$(TYPEDEF)

Defines a list of point vortices.

# Fields

$(TYPEDFIELDS)
"""
struct VortexList
    """list: Array of point vortices."""
    list::Vector{Vortex}
end

VortexList() = VortexList(Vortex[])
VortexList(v...) = VortexList(collect(v))
VortexList(v::VortexList) = v

@propagate_inbounds getindex(A::VortexList, i::Int) = A.list[i]
@propagate_inbounds setindex!(A::VortexList, v::Vortex, i::Int) = A.list[i] = v
@propagate_inbounds getindex(A::VortexList, I...) = A.list[I...]
@propagate_inbounds setindex!(A::VortexList, v, I...) = A.list[I...] = v

iterate(A::VortexList) = iterate(A.list)
iterate(A::VortexList,I) = iterate(A.list,I)
size(A::VortexList) = size(A.list)
length(A::VortexList) = length(A.list)
lastindex(A::VortexList) = lastindex(A.list)

push!(vl::VortexList,v::Vortex) = push!(vl.list,v)

function vcat(vl1::VortexList,vl2::VortexList)
    return VortexList(vcat(vl1.list,vl2.list))
end

"""
$(TYPEDSIGNATURES)

Collect the inertial-space coordinates and strengths of all of the Lagrange points comprising
the vortices in vortex list `vl` and return each assembled set as a vector.
"""
function collect(vl::VortexList)
    xtmp = Float64[]
    ytmp = Float64[]
    Γtmp = Float64[]
    for v in vl
        append!(xtmp,v.x)
        append!(ytmp,v.y)
        append!(Γtmp,v.Γ)
    end
    return xtmp,ytmp,Γtmp
end

collect(vortex::Vortex) = collect(VortexList([vortex]))

"""
$(TYPEDSIGNATURES)

Returns the positions of all point vortices in `vortices` as `VectorData`.
"""
function getpositions(vortices::VortexList)
    positions = VectorData((v->v.x).(vortices.list),(v->v.y).(vortices.list))
    return positions
end

"""
$(TYPEDSIGNATURES)

Returns the strengths of all point vortices in `vortices` as `ScalarData`.
"""
function getstrengths(vortices::VortexList)
    strengths = ScalarData((v->v.Γ).(vortices.list))
    return strengths
end

"""
$(TYPEDSIGNATURES)

Sets the `x` and `y` fields of the point vortices in `vortices` to the entries of `xpositions` and `ypositions`, respectively.
"""
function setpositions!(vortices::VortexList,xpositions,ypositions)
    @assert length(xpositions) == length(vortices)
    @assert length(ypositions) == length(vortices)
    setposition!.(vortices.list,xpositions,ypositions)
end

@recipe f(vl::VortexList) = map(p->p.x,vl.list), map(p->p.y,vl.list)
@recipe f(vl::Vector{Vortex}) = map(p->p.x,vl), map(p->p.y,vl)
