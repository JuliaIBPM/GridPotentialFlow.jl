import Base: *, /

export SuctionParameter, SuctionParameterRange, ModelParameters

const SuctionParameter = Float64

# TODO: use Parameters.jl
struct ModelParameters
    Ub
    U∞
    Γb
    σ
end

function ModelParameters(;Ub=(0.0,0.0), U∞=(0.0,0.0), Γb=nothing, σ=nothing)
    return ModelParameters(Ub,U∞,Γb,σ)
end

struct SuctionParameterRange
    min::SuctionParameter
    max::SuctionParameter

    function SuctionParameterRange(v1,v2)
        if v1 ≤ v2
            min = v1
            max = v2
        else
            min = v2
            max = v1
        end
        new(min,max)
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
