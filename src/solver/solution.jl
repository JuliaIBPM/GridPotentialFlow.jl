abstract type AbstractIBPoissonSolution{TU,TF} end

struct IBPoissonSolution{TU,TF} <: AbstractIBPoissonSolution{TU,TF}
    ψ::TU
    f::TF
end

struct ConstrainedIBPoissonSolution{TU,TF,T} <: AbstractIBPoissonSolution{TU,TF}
    ψ::TU
    f::TF
    ψ₀::Vector{T}
    δΓ_vec::Vector{T}
end

function _scaletophysicalspace!(sol::IBPoissonSolution, Δx::Real)
    sol.ψ .*= Δx # The computed ψ is approximately equal to the continuous streamfunction divided by ∆x.
    sol.f .*= Δx # Because Δψ + Rf = -w, f also has to be scaled back. The result is f = γ*Δs
end

function _scaletophysicalspace!(sol::ConstrainedIBPoissonSolution, Δx::Real)
    sol.ψ .*= Δx # The computed ψ is approximately equal to the continuous streamfunction divided by ∆x.
    sol.f .*= Δx # Because Δψ + Rf = -w, f also has to be scaled back. The result is f = γ*Δs
    sol.ψ₀ .*= Δx
    sol.δΓ_vec .*= Δx
end
