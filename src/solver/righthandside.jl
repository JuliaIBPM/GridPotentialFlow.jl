struct PotentialFlowRHS{T,TU,TF}
    w::TU
    ψb::TF
    f̃limit_kvec::Union{Nothing,Vector{T}}
    f̃min_kvec::Union{Nothing,Vector{T}}
    f̃max_kvec::Union{Nothing,Vector{T}}
end

function PotentialFlowRHS(w::AbstractMatrix,ψb::AbstractVector)
    return PotentialFlowRHS{eltype(w),typeof(w),typeof(ψb)}(w,ψb,nothing,nothing,nothing)
end
