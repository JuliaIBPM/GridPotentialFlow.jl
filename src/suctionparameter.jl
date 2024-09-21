import Base: *, /

export SuctionParameter, SuctionParameterRange, ModelParameters

const SuctionParameter = Float64

const f̃Limit = Float64

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

struct f̃Limits
    min::f̃Limit
    max::f̃Limit

    function f̃Limits(f̃₁,f̃₂)
        f̃sorted = sort([f̃₁,f̃₂])
        new(f̃sorted[1],f̃sorted[2])
    end
end

# Multiply and divide by a constant
function (*)(p::f̃Limits,c::Number)
    return f̃Limits(c*p.min,c*p.max)
end

(*)(c::Number,p::f̃Limits) = *(p,c)

function (/)(p::f̃Limits,c::Number)
    return f̃Limits(p.min/c, p.max/c)
end

function _computef̃limit(SP::SuctionParameter, plate::Polygon, Γ₀)
    f̃ = -SP*2π*platelen(plate)/Γ₀
    return f̃
end

function _computef̃limit(SP::SuctionParameter, body::Body, Γ₀)
    if iszero(SP)
        return 0.0
    else
        throw(ArgumentError("No conversion from suction parameter to limit on f̃ for $(typeof(body)) implemented"))
    end
end

function _computef̃limit(SP::SuctionParameterRange, args...)
    if SP.σmax == SP.σmin
        _computef̃limit(SP.σmax, args...)
    else
        throw(ArgumentError("Ambiguous SuctionParameterRange provided for steady state regularization"))
    end
end

function _computef̃range(SP::SuctionParameter, args...)
    return _computef̃range(SuctionParameterRange(-SP,SP), args...)
end

function _computef̃range(SP::SuctionParameterRange, args...)
    f̃min = _computef̃limit(SP.σmax, args...)
    f̃max = _computef̃limit(SP.σmin, args...)
    return f̃Limits(f̃min,f̃max)
end
