export PotentialFlowRHS

struct UnregularizedPotentialFlowRHS{TU,TF}
    w::TU
    ψb::TF
    Γb::Vector{Float64}
end

struct SteadyRegularizedPotentialFlowRHS{TU,TF,TSP<:Union{SuctionParameter,SuctionParameterRange}}
    w::TU
    ψb::TF
    f̃lim_kvec::Vector{TSP}
end

struct UnsteadyRegularizedPotentialFlowRHS{TU,TF,TSP<:Union{SuctionParameter,SuctionParameterRange}}
    w::TU
    ψb::TF
    f̃lim_kvec::Vector{TSP}
    Γw::Float64
end

function PotentialFlowRHS(w::AbstractMatrix,ψb::AbstractVector;Γ::Union{Nothing,Real,Vector{<:Real}}=nothing)
    if isnothing(Γ)
        Γ = -sum(w)
    end
    return UnregularizedPotentialFlowRHS(w,ψb,[Γ...])
end

function PotentialFlowRHS(w::AbstractMatrix,ψb::AbstractVector,f̃lim_kvec::Vector{TSP}) where {TSP<:Union{SuctionParameter,SuctionParameterRange}}
    return SteadyRegularizedPotentialFlowRHS(w,ψb,f̃lim_kvec)
end

function PotentialFlowRHS(w::AbstractMatrix,ψb::AbstractVector,f̃lim_kvec::Vector{TSP},Γw::Real) where {TSP<:Union{SuctionParameter,SuctionParameterRange}}
    return UnsteadyRegularizedPotentialFlowRHS(w,ψb,f̃lim_kvec,Γw)
end
