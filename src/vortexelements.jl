export PointVortex

struct VortexCache{RCT,ECT} <: ImmersedLayers.AbstractExtraILMCache
    Rcurl :: RCT
    Ecurl :: ECT
end

function VortexCache(base_cache::BasicILMCache)

    # VortexCache(f₀,f̃temp,R̃curl)
end

"""
$(TYPEDEF)

Defines a point vortex with x-position `x`, y-position `y`, and strength `Γ`.

# Fields

$(TYPEDFIELDS)
"""
mutable struct PointVortex
    """x: x-coordinate of the vortex position."""
    x::Float64
    """y: y-coordinate of the vortex position."""
    y::Float64
    """Γ: Strength of the vortex. Positive if counter-clockwise."""
    Γ::Float64
end

"""
$(TYPEDSIGNATURES)

Sets the `x` and `y` fields of the point vortices in the `StructArray` `vl` to the entries of `xnew` and `ynew`, respectively.
"""
function setpositions!(vl::StructArray{PointVortex}, xnew, ynew)
    vl.x .= xnew
    vl.y .= ynew
end

"""
$(TYPEDSIGNATURES)

Sets the `Γ` field of the point vortices in the `StructArray` `vl` to the entries of `Γnew`.
"""
function setstrengths!(vl::StructArray{PointVortex}, Γnew, idx=nothing)
    if isnothing(idx)
        vl.Γ .= Γnew
    else
        setindex!(vl.Γ, Γnew, idx)
    end
end
