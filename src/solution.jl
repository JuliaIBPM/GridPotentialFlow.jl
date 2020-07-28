struct PotentialFlowSolution{T,TU,TF,TB}
    ψ::TU
    f̃::TF
    ψ₀::TB
    δΓ_kvec::Vector{T}
end
