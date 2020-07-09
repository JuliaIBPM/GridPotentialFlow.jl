import LinearAlgebra: I, \, ldiv!

struct PotentialFlowSystem{T,Nk}
    S::SaddleSystem{T}
    f₀::AbstractVector{T}
    ekvec::AbstractVector{<:AbstractVector}
    dkvec::AbstractVector{<:AbstractMatrix}
    f̃kvec::AbstractVector{<:AbstractVector}
end

function PotentialFlowSystem(S::SaddleSystem{T}, f₀::AbstractVector{T}, ekvec::AbstractVector{<:AbstractVector}, dkvec::AbstractVector{<:AbstractMatrix}) where T
    N = length(f₀)
    Nk = length(ekvec)

    # ONLY WORKS FOR A SINGLE EDGE RELEASE FOR NOW!
    @assert Nk < 2 "Only single edge vortex release implemented!"

#     fvec = Array{ScalarData{N,T,Array{T,1}},1}(undef, Nk)
    f̃kvec = fill(ScalarData(N),Nk)

    return PotentialFlowSystem{T,Nk}(S,f₀,ekvec,dkvec,f̃kvec)
end

function ldiv!(sol::Tuple{AbstractMatrix, AbstractVector, AbstractVector, AbstractVector}, sys::PotentialFlowSystem{T,Nk}, rhs::Tuple{AbstractMatrix, AbstractVector, AbstractVector, Real}) where {T,Nk}

    @unpack S, f₀, ekvec, dkvec, f̃kvec = sys
    ψ, f̃, ψ₀, δΓkvec = sol
    negw, ψb, f̃limk, negΓw = rhs

    length(ekvec) == length(dkvec) == length(f̃kvec) == length(δΓkvec) || error("Incompatible number of vortex sheet strength constraints")
    size(negw) == size(ψ) || error("Incompatible number of grid points")

    # Streamfunction and vortex sheet strength without shedded vortices
    tempsol = S\SaddleVector(negw,ψb);
    ψ .= state(tempsol)
    f̃ .= constraint(tempsol)

    N = length(f₀)

    zerovec = ScalarData(N)
    zerovec .= 0.0
    onevec = ScalarData(N)
    onevec .= 1.0

    activef̃limk = fill(0.0,Nk)

    releasevortex = fill(Bool(0),Nk)

    for i in 1:length(ekvec)
        f̃kvec[i] = constraint(S\SaddleVector(-dkvec[i],zerovec));
        -f̃limk[i] ≤ ekvec[i]'*f̃ ≤ f̃limk[i] ? (releasevortex[i] = false; δΓkvec[i] = 0.0) : (releasevortex[i] = true; ekvec[i]'*f̃ < -f̃limk[i] ? activef̃limk[i] = -f̃limk[i] : activef̃limk[i] = f̃limk[i])
    end

    # TODO: IMPLEMENT MULTIPLE EDGES
    # TODO: IMPLEMENT MULTIPLE BODIES
    k = 1
    δΓkvec[k] = -(-negΓw+f₀'*(I - onevec*ekvec[k]')*f̃+sum(f₀)*activef̃limk[k])/(1 + f₀'*(I - onevec*ekvec[k]')*f̃kvec[k]);
    ψ₀[k] = ekvec[k]'*(f̃ + δΓkvec[k]*f̃kvec[k]) - activef̃limk[k]; # is ψ₀
    f̃ .= (I - onevec*ekvec[k]')*(f̃ + δΓkvec[k]*f̃kvec[k]) + onevec*activef̃limk[k];

    ψ .= reshape(-S.A⁻¹*(-reshape(negw,:) + S.B₁ᵀ*f̃ + reshape(dkvec[k]*δΓkvec[k],:)),size(ψ))

    return sol
end

function (\)(sys::PotentialFlowSystem{T,Nk},rhs::Tuple{AbstractMatrix, AbstractVector, AbstractVector, Real}) where {T,Nk}
    negw, ψb, f̃limk, negΓw = rhs
    # ONLY WORKS FOR A SINGLE BODY FOR NOW!
    sol = (similar(negw),similar(ψb),Array{T}(undef, 1),Array{T}(undef, Nk))
    ldiv!(sol,sys,rhs)
    return sol
end
