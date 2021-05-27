export Vortex

"""
$(TYPEDEF)

Defines a point vortex with x-position `x`, y-position `y`, and strength `Γ`.

# Fields

$(TYPEDFIELDS)
"""
mutable struct Vortex
    """x: x-coordinate of the vortex position."""
    x::Float64
    """y: y-coordinate of the vortex position."""
    y::Float64
    """Γ: Strength of the vortex. Positive if counter-clockwise."""
    Γ::Float64
end

"""
$(TYPEDSIGNATURES)

Returns the strengths of all point vortices in the `StructArray` `vl` as `ScalarData`.
"""
function getstrengths(vl::StructArray{Vortex})
    return ScalarData(vl.Γ)
end

"""
$(TYPEDSIGNATURES)

Returns the positions of all point vortices in the `StructArray` `vl` as `VectorData`.
"""
function getpositions(vl::StructArray{Vortex})
    return VectorData(vl.x,vl.y)
end

"""
$(TYPEDSIGNATURES)

Sets the `x` and `y` fields of the point vortices in the `StructArray` `vl` to the entries of `xnew` and `ynew`, respectively.
"""
function setpositions!(vl::StructArray{Vortex}, xnew, ynew)
    vl.x .= xnew
    vl.y .= ynew
end

"""
$(TYPEDSIGNATURES)

Sets the `Γ` field of the point vortices in the `StructArray` `vl` to the entries of `Γnew`.
"""
function setstrengths!(vl::StructArray{Vortex}, Γnew, idx=nothing)
    if isnothing(idx)
        vl.Γ .= Γnew
    else
        vl.Γ[idx] .= Γnew
    end
end
