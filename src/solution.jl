struct PotentialFlowSolution{T,TU,TF}
    ψ::TU
    f̃::TF
    ψ₀::Union{Nothing,Vector{T}}
    δΓ_kvec::Union{Nothing,Vector{T}}
end

function PotentialFlowSolution(ψ::AbstractMatrix,f̃::AbstractVector)
    return PotentialFlowSolution{eltype(ψ),typeof(ψ),typeof(f̃)}(ψ,f̃,nothing,nothing)
end
