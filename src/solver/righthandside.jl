export PotentialFlowRHS

struct UnregularizedPotentialFlowRHS{TU,TF}
    w::TU
    ψb::TF
    Γb::Vector{Float64}
end

struct RegularizedPotentialFlowRHS{TU,TF,TSP<:Union{SuctionParameter,SuctionParameterRange}}
    w::TU
    ψb::TF
    f̃lim_kvec::Vector{TSP}
end

function PotentialFlowRHS(w::AbstractMatrix,ψb::AbstractVector,Γ::Union{Nothing,Float64}=nothing)
    if isnothing(Γ)
        Γ = -sum(w)
    end
    return UnregularizedPotentialFlowRHS(w,ψb,[Γ])
end

function PotentialFlowRHS(w::AbstractMatrix,ψb::AbstractVector,f̃lim_kvec::Vector{TSP}) where {TSP<:Union{SuctionParameter,SuctionParameterRange}}
    return RegularizedPotentialFlowRHS(w,ψb,f̃lim_kvec)
end

function PotentialFlowRHS(w::AbstractMatrix,ψb::AbstractVector,f̃lim_kvec::Tuple{Vector,Vector})
    return RegularizedPotentialFlowRHS(w,ψb,SuctionParameterRange.(f̃lim_kvec[1],f̃lim_kvec[2]))
end

function PotentialFlowRHS(w::AbstractMatrix,ψb::AbstractVector,f̃lim_kvec::Vector{<:Tuple{<:Real,<:Real}})
    return RegularizedPotentialFlowRHS(w,ψb,[SuctionParameterRange(f̃lim[1],f̃lim[2]) for f̃lim in f̃lim_kvec])
end

function PotentialFlowRHS(w::AbstractMatrix,ψb::AbstractVector,f̃lim1_kvec::Vector,f̃lim2_kvec::Vector)
    return RegularizedPotentialFlowRHS(w,ψb,SuctionParameterRange.(f̃lim1_kvec,f̃lim2_kvec))
end
