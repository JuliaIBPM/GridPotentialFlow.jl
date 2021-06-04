abstract type AbstractIBPoissonRHS{TU,TF} end

struct PoissonRHS{TU}
    w::TU
end

struct IBPoissonRHS{TU,TF} <: AbstractIBPoissonRHS{TU,TF}
    w::TU
    ψb::TF
end

struct ConstrainedIBPoissonRHS{T,TU,TF} <: AbstractIBPoissonRHS{TU,TF}
    w::TU
    ψb::TF
    fconstraintRHS::Vector{T}
end

struct UnsteadyRegularizedIBPoissonRHS{TU,TF} <: AbstractIBPoissonRHS{TU,TF}
    w::TU
    ψb::TF
    f̃lim_vec::Vector{f̃Limits}
    Γw_vec::Vector{Float64}
end

function _scaletoindexspace!(rhs::PoissonRHS, Δx::Real)
    rhs.w .*= Δx
end

function _scaletoindexspace!(rhs::IBPoissonRHS, Δx::Real)
    rhs.w .*= Δx
    rhs.ψb ./= Δx
end

function _scaletoindexspace!(rhs::ConstrainedIBPoissonRHS, Δx::Real)
    rhs.w .*= Δx
    rhs.ψb ./= Δx
    rhs.fconstraintRHS ./= Δx
end

function _scaletoindexspace!(rhs::UnsteadyRegularizedIBPoissonRHS, Δx::Real)
    rhs.w .*= Δx
    rhs.ψb ./= Δx
    rhs.f̃lim_vec ./= Δx  # Use same scaling for σ as for f
    rhs.Γw_vec ./= Δx
end
