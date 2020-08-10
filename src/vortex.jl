"""
    Vortex

Vortex type, with x-position `x`, y-position `y`, and strength `Γ`.
"""
mutable struct Vortex
    # Position
    x::Float64
    y::Float64
    # Strength
    Γ::Float64
end

"""
    updateposition!(vortex::Vortex,u::Real,v::Real,Δt::Real)

Update the `x` and `y` fields of `vortex`, given the x-velocity `u`, y-velocity `v`, and timestep `Δt`.
"""
function updateposition!(vortex::Vortex,u::Real,v::Real,Δt::Real)
    vortex.x += Δt*u
    vortex.y += Δt*v
    return vortex
end
