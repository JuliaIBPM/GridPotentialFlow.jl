import Base: *, /

const f̃Limit = Float64

struct ConstraintsCache{VST,VFT,VB22T,ST} <: ImmersedLayers.AbstractExtraILMCache
    ψ₀_vec :: VST
    f₀_vec :: VFT
    B₂₂_vec :: VB22T
    S :: ST
end

function ConstraintsCache(RTLinvR,base_cache::BasicILMCache)
    @unpack bl, gcurl_cache = base_cache
    Nb = length(bl)

    ψ₀_vec = [zeros_gridcurl(base_cache) for i in 1:Nb]
    f₀_vec = [zeros_surfacescalar(base_cache) for i in 1:Nb]
    B₂₂_vec = [zeros_surfacescalar(base_cache) for i in 1:Nb]
    S = zeros(Nb,Nb)

    ones_i = zeros_surfacescalar(base_cache)
    for i in 1:length(bl)
        ones_i .= 0.0
        ones_i[getrange(bl,i)] .= 1.0
        f₀_vec[i] .= -RTLinvR\ones_i  #Should use ldiv!(γtemp,RTLinvR) for no allocation or even better implement CG so we don't have to construct RTLinvR.
        regularize!(gcurl_cache,f₀_vec[i],base_cache)
        inverse_laplacian!(ψ₀_vec[i],gcurl_cache,base_cache)
        ψ₀_vec[i] .*= -1
    end

    # Fill up B₂₂ matrix. Not correct for shedding bodies
    for i in 1:length(bl)
        B₂₂_vec[i][getrange(bl,i)] .= createB₂₂entry(bl[i],f₀_vec[i][getrange(bl,i)],base_cache)
    end

    # Precompute S matrix. Not correct for shedding bodies
    for i in 1:Nb, j in 1:Nb
        S[i,j] = -B₂₂_vec[i]'*f₀_vec[j]
    end

    ConstraintsCache{typeof(ψ₀_vec),typeof(f₀_vec),typeof(B₂₂_vec),typeof(S)}(ψ₀_vec,f₀_vec,B₂₂_vec,S)
end

function createB₂₂entry(body::PotentialFlowBody{N,C,ST,0},f₀,base_cache) where {N,C,ST}
    B₂₂entry = ScalarData(N)
    B₂₂entry .= 1.0
    return B₂₂entry
end

function createB₂₂entry(body::PotentialFlowBody{N,C,ST,1},f₀,base_cache) where {N,C,ST}
    B₂₂entry = ScalarData(N)
    B₂₂entry[body.edges[1]] = 1.0
    B₂₂entry ./= f₀
    return B₂₂entry
end

function createB₂₂entry(body::Body,f₀,base_cache)
    return createB₂₂entry(PotentialFlowBody(body),f₀,base_cache)
end

function constraint(body::PotentialFlowBody{N,C,ST,0}) where {N,C,ST}
    return body.Γ
end

function constraint(body::PotentialFlowBody{N,C,ST,1}) where {N,C,ST}
    return _computef̃limit(body.σ[1],body,nothing)
end

function constraint(body::Body) where {N,C,ST}
    return 0.0
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

# function _computef̃limit(SP::SuctionParameter, plate::Plate, Γ₀)
#     f̃ = -SP*2π*plate.len/Γ₀
#     return f̃
# end

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
