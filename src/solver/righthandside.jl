export PotentialFlowRHS


struct BasicPotentialFlowRHS{TU,TF}
    w::TU
    ψb::TF
end









const f̃Limit = Float64

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

struct UnregularizedPotentialFlowRHS{TU,TF}
    w::TU
    ψb::TF
    Γb::Vector{Float64}
end

struct SteadyRegularizedPotentialFlowRHS{TU,TF}
    w::TU
    ψb::TF
    f̃lim_kvec::Vector{f̃Limit}
end

mutable struct UnsteadyRegularizedPotentialFlowRHS{TU,TF}
    w::TU
    ψb::TF
    f̃lim_kvec::Vector{f̃Limits}
    Γw::Float64
end

function PotentialFlowRHS(w::AbstractMatrix,ψb::AbstractVector;Γ::Union{Nothing,Real,Vector{<:Real}}=nothing)
    if isnothing(Γ)
        Γ = -sum(w)
    end
    return UnregularizedPotentialFlowRHS(w,ψb,[Γ...])
end

function PotentialFlowRHS(w::AbstractMatrix,ψb::AbstractVector,f̃lim_kvec::Vector{f̃Limit})
    return SteadyRegularizedPotentialFlowRHS(w,ψb,f̃lim_kvec)
end

function PotentialFlowRHS(w::AbstractMatrix,ψb::AbstractVector,f̃lim_kvec::Vector{f̃Limits},Γw::Real)
    return UnsteadyRegularizedPotentialFlowRHS(w,ψb,f̃lim_kvec,Γw)
end

function _scaletoindexspace!(rhs::UnregularizedPotentialFlowRHS,Δx::Real)
    rhs.w .*= Δx
    rhs.ψb ./= Δx
    if !isnothing(rhs.Γb)
        rhs.Γb ./= Δx
    end
end

function _scaletoindexspace!(rhs::SteadyRegularizedPotentialFlowRHS,Δx::Real)
    rhs.w .*= Δx
    rhs.ψb ./= Δx
    rhs.f̃lim_kvec ./= Δx
end

function _scaletoindexspace!(rhs::UnsteadyRegularizedPotentialFlowRHS,Δx::Real)
    rhs.w .*= Δx
    rhs.ψb ./= Δx
    rhs.f̃lim_kvec ./= Δx  # Use same scaling for σ as for f
    rhs.Γw /= Δx
end
