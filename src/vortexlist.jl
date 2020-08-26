import Base: @propagate_inbounds,getindex, setindex!,iterate,size,length,push!,
              collect, view, vcat, lastindex

"""
    VortexList

"""
struct VortexList
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

push!(bl::VortexList,b::Vortex) = push!(bl.list,b)

function vcat(vl1::VortexList,vl2::VortexList)
    return VortexList(vcat(vl1.list,vl2.list))
end

"""
    collect(vl::vortexlist) -> Vector{Float64}, Vector{Float64}
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


function getpositions(vortices::VortexList)
    positions = VectorData((v->v.x).(vortices.list),(v->v.y).(vortices.list))
    return positions
end

function getstrengths(vortices::VortexList)
    strengths = ScalarData((v->v.Γ).(vortices.list))
    return strengths
end

function setpositions!(vortices::VortexList,xpositions,ypositions)
    setposition!.(vortices.list,xpositions,ypositions)
end
