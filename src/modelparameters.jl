import Base: *, /

export SuctionParameter, SuctionParameterRange, ModelParameters

const SuctionParameter = Float64

struct SuctionParameterRange
    σmin::SuctionParameter
    σmax::SuctionParameter

    function SuctionParameterRange(σ₁,σ₂)
        σsorted = sort([σ₁,σ₂])
        new(σsorted[1],σsorted[2])
    end
end

# Multiply and divide by a constant
function (*)(p::SuctionParameterRange,c::Number)
  return SuctionParameterRange(c*p.min,c*p.max)
end

(*)(c::Number,p::SuctionParameterRange) = *(p,c)

function (/)(p::SuctionParameterRange,c::Number)
  return SuctionParameterRange(p.min/c, p.max/c)
end

"""
$(TYPEDEF)

Defines the model parameters for `VortexModel`.

# Fields

$(TYPEDFIELDS)
"""
struct ModelParameters
    """Ub: Array that contains the velocity of each body or `nothing` if there is no body. Each entry should be a collection of the x- and y-components of the velocity."""
    Ub::Union{AbstractArray,Nothing}
    """U∞: Collection of the x- and y-components of the free stream velocity."""
    U∞
    """Γb: Array that contains the desired circulation of each body or `nothing` if there is no body. This will only be enforced if the body has no regularized edges."""
    Γb::Union{AbstractArray,Nothing}
    """σ: Array that contains the suction parameter or suction parameter range for each edge or `nothing` if there are no regularized edges or if the Kutta condition is to be enforced."""
    σ::Union{AbstractArray{Union{SuctionParameter,SuctionParameterRange}},Nothing}
end

"""
$(TYPEDSIGNATURES)

Construct the parameters for a vortex model using the given function.
"""
function ModelParameters(;Ub=nothing, U∞=(0.0,0.0), Γb=nothing, σ=nothing)
    return ModelParameters(Ub,U∞,Γb,σ)
end
