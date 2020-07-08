"""
    Vortex{T}

Vortex type, with x-position `x`, y-position `y`, and strength `Γ`, all of type `T`.
"""
mutable struct Vortex{T}
    # Position
    x::T
    y::T
    # Strength
    Γ::T
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
