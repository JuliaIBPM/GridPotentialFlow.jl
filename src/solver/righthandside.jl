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

mutable struct UnsteadyRegularizedPotentialFlowRHS{TU,TF,TSP<:Union{SuctionParameter,SuctionParameterRange}}
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

function _scaletoindexspace!(rhs::UnregularizedPotentialFlowRHS,Δx::Real)
    rhs.w .*= Δx
    rhs.ψb ./= Δx
    if !isnothing(rhs.Γb)
        rhs.Γb ./= Δx
    end
end

function _scaletoindexspace!(rhs::SteadyRegularizedPotentialFlowRHS,Δx::Real)
    rhs.w .*= Δx
    rhs.ψb ./= Δx
    rhs.f̃lim_kvec ./= Δx
end

function _scaletoindexspace!(rhs::UnsteadyRegularizedPotentialFlowRHS,Δx::Real)
    rhs.w .*= Δx
    rhs.ψb ./= Δx
    rhs.f̃lim_kvec ./= Δx  # Use same scaling for σ as for f
    rhs.Γw /= Δx
end
