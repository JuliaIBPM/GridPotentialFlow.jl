import LinearAlgebra: I, \, ldiv!

# abstract type PotentialFlowSystem end
# abstract type SystemType end
#
# abstract type NonRegularized <: SystemType end
# abstract type SteadyRegularized <: SystemType end
# abstract type UnsteadyRegularized <: SystemType end

struct UnregularizedPotentialFlowSystem
    S::SaddleSystem
end

struct RegularizedPotentialFlowSystem{TU,TF}
    S::SaddleSystem
    f₀::TF
    ek_vec::Vector{BodyUnitVector{TF}}
    dk_vec::Union{Nothing,Vector{TU}}
    f̃k_vec::Union{Nothing,Vector{TF}}
end

function PotentialFlowSystem(S::SaddleSystem{T,Ns,Nc,TU,TF}) where {T,Ns,Nc,TU,TF}
    return UnregularizedPotentialFlowSystem(S)
end

function PotentialFlowSystem(S::SaddleSystem{T,Ns,Nc,TU,TF}, f₀::TF, ek_vec::Vector{BodyUnitVector{TF}}) where {T,Ns,Nc,TU,TF}
    return RegularizedPotentialFlowSystem{TU,TF}(S,f₀,ek_vec,nothing,nothing)
end

function PotentialFlowSystem(S::SaddleSystem{T,Ns,Nc,TU,TF}, f₀::TF, ek_vec::Vector{BodyUnitVector{TF}}, dk_vec::Vector{TU}) where {T,Ns,Nc,TU,TF}
    size(ek_vec) == size(dk_vec) || error("Incompatible number of elements in ekvec and dkvec")
    Nk = length(ek_vec)
    f̃k_vec = Vector{TF}(undef, Nk)
    return RegularizedPotentialFlowSystem{TU,TF}(S, f₀, ek_vec, dk_vec, f̃k_vec)
end

function ldiv!(sol::PotentialFlowSolution{TU,TF}, sys::RegularizedPotentialFlowSystem{TU,TF}, rhs::PotentialFlowRHS{TU,TF}) where {TU,TF}

    # TODO: maybe type check here instead of in arguments list?

    # # Streamfunction and vortex sheet strength without shedded vortices
    # tempsol = S\SaddleVector(negw,ψb);
    # ψ .= state(tempsol)
    # f̃ .= constraint(tempsol)
    #
    # N = length(f₀)
    #
    # zerovec = ScalarData(N)
    # zerovec .= 0.0
    # onevec = ScalarData(N)
    # onevec .= 1.0
    #
    # activef̃limk = fill(0.0,Nk)
    #
    # releasevortex = fill(Bool(0),Nk)
    #
    # for i in 1:length(ekvec)
    #     f̃kvec[i] = constraint(S\SaddleVector(-dkvec[i],zerovec));
    #     -f̃limk[i] ≤ ekvec[i]'*f̃ ≤ f̃limk[i] ? (releasevortex[i] = false; δΓkvec[i] = 0.0) : (releasevortex[i] = true; ekvec[i]'*f̃ < -f̃limk[i] ? activef̃limk[i] = -f̃limk[i] : activef̃limk[i] = f̃limk[i])
    # end
    #
    # # TODO: IMPLEMENT MULTIPLE EDGES
    # # TODO: IMPLEMENT MULTIPLE BODIES
    # k = 1
    # δΓkvec[k] = -(-negΓw+f₀'*(I - onevec*ekvec[k]')*f̃+sum(f₀)*activef̃limk[k])/(1 + f₀'*(I - onevec*ekvec[k]')*f̃kvec[k]);
    # ψ₀[k] = ekvec[k]'*(f̃ + δΓkvec[k]*f̃kvec[k]) - activef̃limk[k]; # is ψ₀
    # f̃ .= (I - onevec*ekvec[k]')*(f̃ + δΓkvec[k]*f̃kvec[k]) + onevec*activef̃limk[k];
    #
    # ψ .= reshape(-S.A⁻¹*(-reshape(negw,:) + S.B₁ᵀ*f̃ + reshape(dkvec[k]*δΓkvec[k],:)),size(ψ))

    return sol
end

# function (\)(sys::PotentialFlowSystem{Nk,T},rhs::Tuple{AbstractMatrix, AbstractVector, AbstractVector, Real}) where {Nk,T}
#     negw, ψb, f̃limk, negΓw = rhs
#     # ONLY WORKS FOR A SINGLE BODY FOR NOW!
#     sol = (similar(negw),similar(ψb),Array{T}(undef, 1),Array{T}(undef, Nk))
#     ldiv!(sol,sys,rhs)
#     return sol
# end
