export PotentialFlowRHS

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

struct IBPoissonRHS{TU,TF}
    w::TU
    ψb::TF
end

struct ConstrainedIBPoissonRHS{T,TU,TF}
    w::TU
    ψb::TF
    fconstraintRHS::Vector{T}
end

struct UnsteadyRegularizedIBPoissonRHS{TU,TF}
    w::TU
    ψb::TF
    f̃lim_vec::Vector{f̃Limits}
    Γw_vec::Vector{Float64}
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
