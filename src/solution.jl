struct PotentialFlowSolution{TU,TF,T}
    ψ::TU
    f::TF
    ψ₀::T
    δΓvec::Vector{T}
end
