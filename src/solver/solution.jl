export PotentialFlowSolution

abstract type PotentialFlowSolution end

struct UnregularizedPotentialFlowSolution{T,TU,TF} <: PotentialFlowSolution
    ψ::TU
    f::TF
    ψ₀::Vector{T}
end

struct SteadyRegularizedPotentialFlowSolution{T,TU,TF} <: PotentialFlowSolution
    ψ::TU
    f̃::TF
    ψ₀::Vector{T}
end

struct UnsteadyRegularizedPotentialFlowSolution{T,TU,TF} <: PotentialFlowSolution
    ψ::TU
    f̃::TF
    ψ₀::Vector{T}
    δΓ_kvec::Vector{T}
end

function PotentialFlowSolution(ψ::TU,f::TF) where {TU,TF}
    return UnregularizedPotentialFlowSolution{Float64,TU,TF}(ψ,f,Float64[])
end

function PotentialFlowSolution(ψ::TU,f̃::TF,ψ₀::Vector{T}) where {T,TU,TF}
    return SteadyRegularizedPotentialFlowSolution{T,TU,TF}(ψ,f̃,ψ₀)
end

function PotentialFlowSolution(ψ::TU,f̃::TF,ψ₀::Vector{T},δΓ_kvec::Vector{T}) where {T,TU,TF}
    return UnsteadyRegularizedPotentialFlowSolution{T,TU,TF}(ψ,f̃,ψ₀,δΓ_kvec)
end
