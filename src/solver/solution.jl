abstract type FlowType end

struct UnregularizedPotentialFlowSolution{TU,TF}
    ψ::TU
    f::TF
end

struct SteadyRegularizedPotentialFlowSolution{T,TU,TF}
    ψ::TU
    f̃::TF
    ψ₀::Vector{T}
end

struct UnsteadyRegularizedPotentialFlowSolution{T,TU,TF}
    ψ::TU
    f̃::TF
    ψ₀::Vector{T}
    δΓ_kvec::Union{Nothing,Vector{T}}
end
