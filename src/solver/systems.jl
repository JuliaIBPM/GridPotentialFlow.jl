import LinearAlgebra: I, \, ldiv!
import SparseArrays: AbstractSparseMatrix, SparseMatrixCSC, sparse
import Statistics: mean

abstract type PotentialFlowSystem end

struct UnregularizedPotentialFlowSystem <: PotentialFlowSystem
    S::SaddleSystem
end

struct RegularizedPotentialFlowSystem{Nb,Nk,T,TU,TF,TE} <: PotentialFlowSystem
    S̃::SaddleSystem
    f₀::TF
    e_kvec::Vector{TE}
    d_kvec::Vector{TU}
    f̃_kvec::Vector{TF}
    P_kvec::Vector{AbstractSparseMatrix}
    zeros::TF
    ones::TF
    _w_buf::TU

    function RegularizedPotentialFlowSystem(S̃::SaddleSystem{T,Ns,Nc,TU,TF}, f₀::TF, e_kvec::Vector{TE}, d_kvec::Vector{TU}, f̃_kvec::Vector{TF}) where {T,Ns,Nc,TU,TF,TE}
        _w_buf = TU()
        P_kvec = _computesparsekuttaoperator.(e_kvec)
        zeros = TF()
        zeros .= 0
        ones = TF()
        ones .= 1
        new{1,length(e_kvec),T,TU,TF,TE}(S̃, f₀, e_kvec, d_kvec, f̃_kvec, P_kvec, zeros, ones, _w_buf)
    end

    # TODO check cost benefit of calling this constructor instead of the previous one because this one is unsafe
    function RegularizedPotentialFlowSystem(S̃::SaddleSystem{T,Ns,Nc,TU,TF}, f₀::TF, e_kvec::Vector{TE}, d_kvec::Vector{TU}, f̃_kvec::Vector{TF}, P_kvec::Vector{AbstractSparseMatrix}, zeros::TF, ones::TF, _w_buf::TU) where {T,Ns,Nc,TU,TF,TE}
        new{1,length(e_kvec),T,TU,TF,TE}(S̃, f₀, e_kvec, d_kvec, f̃_kvec, P_kvec, zeros, ones, _w_buf)
    end

end

function PotentialFlowSystem(S::SaddleSystem{T,Ns,Nc,TU,TF}) where {T,Ns,Nc,TU,TF}
    return UnregularizedPotentialFlowSystem(S)
end

function PotentialFlowSystem(S̃::SaddleSystem{T,Ns,Nc,TU,TF}, f₀::TF, e_kvec::Vector{TE}) where {T,Ns,Nc,TU,TF,TE}
    return RegularizedPotentialFlowSystem(S̃,f₀,e_kvec,TU[],TF[])
end

function PotentialFlowSystem(S̃::SaddleSystem{T,Ns,Nc,TU,TF}, f₀::TF, e_kvec::Vector{TE}, d_kvec::Vector{TU}) where {T,Ns,Nc,TU,TF,TE}
    f̃_kvec = fill(TF(), length(e_kvec))
    return RegularizedPotentialFlowSystem(S̃,f₀,e_kvec,d_kvec,f̃_kvec)
end

function ldiv!(sol::UnregularizedPotentialFlowSolution{TU,TF}, sys::UnregularizedPotentialFlowSystem, rhs::UnregularizedPotentialFlowRHS{TU,TF}) where {T,TU,TF,TE}

    @unpack S = sys
    @unpack ψ, f = sol
    @unpack w, ψb = rhs

    temp_sol = S\SaddleVector(-w,ψb)
    ψ .= state(temp_sol)
    f .= constraint(temp_sol)

    return sol
end

function ldiv!(sol::SteadyRegularizedPotentialFlowSolution{T,TU,TF}, sys::RegularizedPotentialFlowSystem{Nb,Nk,T,TU,TF,TE}, rhs::RegularizedPotentialFlowRHS{TU,TF,SuctionParameter}) where {Nb,Nk,T,TU,TF,TE,TSP}

    @unpack S̃, f₀, e_kvec, P_kvec, zeros, ones = sys
    @unpack ψ, f̃, ψ₀ = sol
    @unpack w, ψb, f̃lim_kvec = rhs

    # TODO: IMPLEMENT MULTIPLE BODIES

    @assert size(f̃lim_kvec) == size(e_kvec)

    f̃ .= constraint(S̃\SaddleVector(-w,ψb))
    ψ₀ .= mean(e_kvec)'*f̃ .- mean(f̃lim_kvec)
    f̃ .= mean(P_kvec)*f̃ + ones*mean(f̃lim_kvec)
    ψ .= reshape(-S̃.A⁻¹*(reshape(w,:) + S̃.B₁ᵀ*f̃),size(w))

    return sol
end

function ldiv!(sol::UnsteadyRegularizedPotentialFlowSolution{T,TU,TF}, sys::RegularizedPotentialFlowSystem{Nb,Nk,T,TU,TF,TE}, rhs::RegularizedPotentialFlowRHS{TU,TF,TSP}) where {Nb,Nk,T,TU,TF,TE,TSP}

    @unpack S̃, f₀, e_kvec, d_kvec, f̃_kvec, P_kvec, zeros, ones, _w_buf = sys
    @unpack ψ, f̃, ψ₀, δΓ_kvec = sol
    @unpack w, ψb, f̃lim_kvec = rhs

    # TODO: IMPLEMENT MULTIPLE BODIES
    @assert length(d_kvec) == Nk
    @assert length(f̃_kvec) == Nk
    @assert length(f̃lim_kvec) == Nk
    @assert length(δΓ_kvec) == Nk

    for k in 1:Nk
        f̃_kvec[k] = constraint(S̃\SaddleVector(-d_kvec[k],zeros)) # Maybe move this to system constructor?
    end

    f̃ .= constraint(S̃\SaddleVector(-w,ψb))
    k_sheddingedges, activef̃lim_kvec = _findsheddingedges(Nk, e_kvec, f̃, f̃lim_kvec)

    # TODO: add loop below to correct vortex strengths (after adding constraint function to ConstrainedSystems)
    # while true
    #     f̃ .= constraint(S̃\SaddleVector(-w .+ sum(δΓ_kvec.*d_kvec,ψb))
    #     k_sheddingedges, activef̃limits_kvec = _findsheddingedges(Nk, e_kvec, f̃, f̃min_kvec, f̃max_kvec)
    #     if k_sheddingedges
    #         break
    #     end
    #     δΓ_kvec .= _computepointvortexstrengths(Nk, k_sheddingedges::Vector{Integer}, P_kvec, f̃_kvec, f̃lim_kvec, f₀, w)
    # end

    δΓ_kvec .= _computepointvortexstrengths(Nk, k_sheddingedges::Vector{Integer}, P_kvec, f̃_kvec, activef̃lim_kvec, f̃, f₀, w)

    _w_buf .= w .+ sum(δΓ_kvec.*d_kvec)

    steadyrhs = PotentialFlowRHS(_w_buf, ψb, activef̃lim_kvec)
    steadysys = RegularizedPotentialFlowSystem(S̃, f₀, e_kvec[k_sheddingedges], TU[], TF[], P_kvec[k_sheddingedges], zeros, ones, _w_buf)
    steadysol = SteadyRegularizedPotentialFlowSolution(ψ,f̃,ψ₀)
    ldiv!(steadysol,steadysys,steadyrhs)

    return sol
end

function (\)(sys::UnregularizedPotentialFlowSystem, rhs::UnregularizedPotentialFlowRHS{TU,TF}) where {TU,TF}
    sol = UnregularizedPotentialFlowSolution(TU(),TF())
    ldiv!(sol,sys,rhs)
    return sol
end

function (\)(sys::RegularizedPotentialFlowSystem{Nb,Nk,T,TU,TF,TE}, rhs::RegularizedPotentialFlowRHS{TU,TF,TSP}) where {Nb,Nk,T,TU,TF,TE,TSP}
    if isempty(sys.d_kvec)
        println("d_kvec not set in system. Providing steady regularized solution.")
        sol = SteadyRegularizedPotentialFlowSolution(TU(),TF(),zeros(T,Nb))
    else
        sol = UnsteadyRegularizedPotentialFlowSolution(TU(),TF(),zeros(T,Nb),zeros(T,Nk))
    end
    ldiv!(sol,sys,rhs)
    return sol
end

function _computesparsekuttaoperator(e::AbstractVector)::SparseMatrixCSC
    return sparse(I - ones(length(e))*e')
end

function _findsheddingedges(Nk, e_kvec, f̃, f̃lim_kvec::Vector{SuctionParameter})

    f̃lim_kvec_range = SuctionParameterRange.(-f̃lim_kvec,f̃lim_kvec)

    return _findsheddingedges(Nk, e_kvec, f̃, f̃lim_kvec_range)

end

function _findsheddingedges(Nk, e_kvec, f̃, f̃lim_kvec::Vector{SuctionParameterRange})

    k_sheddingedges = Vector{Integer}()
    activef̃lim_kvec = Vector{SuctionParameter}()

    for k in 1:Nk
        if e_kvec[k]'*f̃ < f̃lim_kvec[k].min
            push!(activef̃lim_kvec,f̃lim_kvec[k].min)
            push!(k_sheddingedges,k)
        elseif e_kvec[k]'*f̃ > f̃lim_kvec[k].max
            push!(activef̃lim_kvec,f̃lim_kvec[k].max)
            push!(k_sheddingedges,k)
        else
            push!(activef̃lim_kvec,Inf)
        end
    end

    return k_sheddingedges, activef̃lim_kvec

end

function _computepointvortexstrengths(Nk, k_sheddingedges::Vector{Integer}, P_kvec, f̃_kvec, f̃lim_kvec, f̃, f₀, w)

    δΓsys = zeros(eltype(w),Nk,Nk)
    δΓrhs = zeros(eltype(w),Nk)
    Γ₀ = sum(f₀)
    Γw = sum(w)

    for k in 1:Nk
        if k in k_sheddingedges
            for l in k_sheddingedges
                δΓsys[l,k] = 1 + f₀'*P_kvec[l]*f̃_kvec[k]
            end
            δΓrhs[k] = Γw + f₀'*P_kvec[k]*f̃ + Γ₀*f̃lim_kvec[k]
        else
            δΓsys[k,k] = 1
            δΓrhs[k] = 0
        end
    end

    δΓ_kvec = -δΓsys\δΓrhs

    return δΓ_kvec
end

function setd_kvec!(sys::RegularizedPotentialFlowSystem{Nb,Nk,T,TU,TF,TE}, d_kvec::Vector{TU}) where {Nb,Nk,T,TU,TF,TE}

    sys.d_kvec = d_kvec
    sys.f_kvec = fill(TF(), length(d_kvec))

end
