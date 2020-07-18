struct PotentialFlowRHS{TU,TF,T}
    w::TU
    ψb::TF
    f̃min::Vector{T}
    f̃max::Vector{T}
end
