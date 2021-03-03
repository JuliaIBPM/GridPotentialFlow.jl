export PotentialFlowSolution

abstract type PotentialFlowSolution end

struct UnregularizedPotentialFlowSolution{T,TU,TF} <: PotentialFlowSolution
    ψ::TU
    f::TF
    ψ₀::Vector{T}
end

struct SteadyRegularizedPotentialFlowSolution{T,TU,TF} <: PotentialFlowSolution
    ψ::TU
    f::TF
    ψ₀::Vector{T}
end

struct UnsteadyRegularizedPotentialFlowSolution{T,TU,TF} <: PotentialFlowSolution
    ψ::TU
    f::TF
    ψ₀::Vector{T}
    δΓ_kvec::Vector{T}
end

function PotentialFlowSolution(ψ::TU,f::TF) where {TU,TF}
    return UnregularizedPotentialFlowSolution{Float64,TU,TF}(ψ,f,Float64[])
end

function PotentialFlowSolution(ψ::TU,f::TF,ψ₀::Vector{T}) where {T,TU,TF}
    return SteadyRegularizedPotentialFlowSolution{T,TU,TF}(ψ,f,ψ₀)
end

function PotentialFlowSolution(ψ::TU,f::TF,ψ₀::Vector{T},δΓ_kvec::Vector{T}) where {T,TU,TF}
    return UnsteadyRegularizedPotentialFlowSolution{T,TU,TF}(ψ,f,ψ₀,δΓ_kvec)
end

function _scaletophysicalspace!(sol::UnregularizedPotentialFlowSolution,Δx::Real)
    sol.ψ .*= Δx # The computed ψ is approximately equal to the continuous streamfunction divided by ∆x.
    sol.f .*= Δx # Because Δψ + Rf = -w, f also has to be scaled back. The result is f = γ*Δs
    sol.ψ₀ .*= Δx
end

function _scaletophysicalspace!(sol::SteadyRegularizedPotentialFlowSolution,Δx::Real)
    sol.ψ .*= Δx # The computed ψ is approximately equal to the continuous streamfunction divided by ∆x.
    sol.f .*= Δx # Because Δψ + Rf = -w, f also has to be scaled back. The result is f = γ*Δs
    sol.ψ₀ .*= Δx
end
function _scaletophysicalspace!(sol::UnsteadyRegularizedPotentialFlowSolution,Δx::Real)
    sol.ψ .*= Δx # The computed ψ is approximately equal to the continuous streamfunction divided by ∆x.
    sol.f .*= Δx # Because Δψ + Rf = -w, f also has to be scaled back. The result is f = γ*Δs
    sol.ψ₀ .*= Δx
    sol.δΓ_kvec .*= Δx
end
