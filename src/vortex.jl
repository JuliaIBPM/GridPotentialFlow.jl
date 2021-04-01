export Vortex, updateposition!

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
    """Γ: strength of the vortex. Positive if counter-clockwise."""
    Γ::Float64
end

"""
$(TYPEDSIGNATURES)

Updates the `x` and `y` fields of `vortex`, given the x-velocity `u`, y-velocity `v`, and timestep `Δt`.
"""
function updateposition!(vortex::Vortex,u::Real,v::Real,Δt::Real)
    vortex.x += Δt*u
    vortex.y += Δt*v
    return vortex
end

"""
$(TYPEDSIGNATURES)

Sets the `x` and `y` fields of `vortex` to the `xnew` and `ynew`, respectively.
"""
function setposition!(vortex::Vortex,xnew::Real,ynew::Real)
    vortex.x = xnew
    vortex.y = ynew
    return vortex
end
