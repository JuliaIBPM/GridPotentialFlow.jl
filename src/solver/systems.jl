import LinearAlgebra: I, \, ldiv!
import SparseArrays: AbstractSparseMatrix, SparseMatrixCSC, sparse
import Statistics: mean

abstract type FlowType end
abstract type Steady <: FlowType end
abstract type Unsteady <: FlowType end

struct UnregularizedPotentialFlowSystem
    S::SaddleSystem
end

struct RegularizedPotentialFlowSystem{FlowType,Nb,Nk,T,TU,TF,TE}
    S̃::SaddleSystem
    f₀::TF
    e_kvec::Vector{TE}
    d_kvec::Union{Nothing,Vector{TU}}
    f̃_kvec::Union{Nothing,Vector{TF}}
    P_kvec::Vector{AbstractSparseMatrix}
    zeros::TF
    ones::TF
    _w_buf::Union{Nothing,TU}

    function RegularizedPotentialFlowSystem(flowtype::Type{<:FlowType},S̃::SaddleSystem{T,Ns,Nc,TU,TF}, f₀::TF, e_kvec::Vector{TE}, d_kvec::Union{Nothing,Vector{TU}}=nothing, f̃_kvec::Union{Nothing,Vector{TF}}=nothing) where {T,Ns,Nc,TU,TF,TE}
        if flowtype == Unsteady
            if isnothing(d_kvec) error("d_kvec is not defined") end
            if isnothing(f̃_kvec) error("f̃_kvec is not defined") end
            if size(e_kvec) != size(d_kvec) error("Incompatible number of elements in e_kvec and d_kvec") end
            if size(e_kvec) != size(f̃_kvec) error("Incompatible number of elements in e_kvec and f̃_kvec") end
            _w_buf = TU()
        else
            _w_buf = nothing
        end
        P_kvec = _computesparsekuttaoperator.(e_kvec)
        zeros = TF()
        zeros .= 0
        ones = TF()
        ones .= 1
        new{flowtype,1,length(e_kvec),T,TU,TF,TE}(S̃, f₀, e_kvec, d_kvec, f̃_kvec, P_kvec, zeros, ones, _w_buf)
    end

    # TODO check cost benefit of calling this constructor instead of the previous one because this one is unsafe
    function RegularizedPotentialFlowSystem(flowtype::Type{<:FlowType},S̃::SaddleSystem{T,Ns,Nc,TU,TF}, f₀::TF, e_kvec::Vector{TE}, d_kvec::Union{Nothing,Vector{TU}}, f̃_kvec::Union{Nothing,Vector{TF}}, P_kvec::Vector{AbstractSparseMatrix}, zeros::TF, ones::TF, _w_buf::Union{Nothing,TU}) where {T,Ns,Nc,TU,TF,TE}
        new{flowtype,1,length(e_kvec),T,TU,TF,TE}(S̃, f₀, e_kvec, d_kvec, f̃_kvec, P_kvec, zeros, ones, _w_buf)
    end

end

function PotentialFlowSystem(S::SaddleSystem{T,Ns,Nc,TU,TF}) where {T,Ns,Nc,TU,TF}
    return UnregularizedPotentialFlowSystem(S)
end

function PotentialFlowSystem(S̃::SaddleSystem{T,Ns,Nc,TU,TF}, f₀::TF, e_kvec::Vector{TE}) where {T,Ns,Nc,TU,TF,TE}
    return RegularizedPotentialFlowSystem(Steady,S̃,f₀,e_kvec)
end

function PotentialFlowSystem(S̃::SaddleSystem{T,Ns,Nc,TU,TF}, f₀::TF, e_kvec::Vector{TE}, d_kvec::Vector{TU}) where {T,Ns,Nc,TU,TF,TE}
    f̃_kvec = fill(TF(), length(e_kvec))
    return RegularizedPotentialFlowSystem(Unsteady,S̃,f₀,e_kvec,d_kvec,f̃_kvec)
end

function ldiv!(sol::PotentialFlowSolution{T,TU,TF}, sys::UnregularizedPotentialFlowSystem, rhs::PotentialFlowRHS{T,TU,TF}) where {FlowType,T,TU,TF,TE}

    @unpack S = sys
    @unpack ψ, f̃ = sol #TODO: create sol with f instead of f̃
    @unpack w, ψb = rhs

    temp_sol = S\SaddleVector(-w,ψb)
    ψ .= state(temp_sol)
    f̃ .= constraint(temp_sol)

    return sol
end

function ldiv!(sol::PotentialFlowSolution{T,TU,TF}, sys::RegularizedPotentialFlowSystem{Steady,Nb,Nk,T,TU,TF,TE}, rhs::PotentialFlowRHS{T,TU,TF}) where {FlowType,Nb,Nk,T,TU,TF,TE}

    @unpack S̃, f₀, e_kvec, P_kvec, zeros, ones = sys
    @unpack ψ, f̃, ψ₀ = sol
    @unpack w, ψb, f̃limit_kvec = rhs

    # TODO: IMPLEMENT MULTIPLE BODIES

    @assert size(f̃limit_kvec) == size(e_kvec)

    f̃ .= constraint(S̃\SaddleVector(-w,ψb))
    ψ₀ .= mean(e_kvec)'*f̃ .- mean(f̃limit_kvec)
    f̃ .= mean(P_kvec)*f̃ + ones*mean(f̃limit_kvec)
    ψ .= reshape(-S̃.A⁻¹*(reshape(w,:) + S̃.B₁ᵀ*f̃),size(w))

    return sol
end

function ldiv!(sol::PotentialFlowSolution{T,TU,TF}, sys::RegularizedPotentialFlowSystem{Unsteady,Nb,Nk,T,TU,TF,TE}, rhs::PotentialFlowRHS{T,TU,TF}) where {FlowType,Nb,Nk,T,TU,TF,TE}

    @unpack S̃, f₀, e_kvec, d_kvec, f̃_kvec, P_kvec, zeros, ones, _w_buf = sys
    @unpack ψ, f̃, ψ₀, δΓ_kvec = sol
    @unpack w, ψb, f̃min_kvec, f̃max_kvec = rhs

    # TODO: IMPLEMENT MULTIPLE BODIES
    @assert f̃min_kvec != nothing
    @assert f̃max_kvec != nothing
    @assert length(f̃min_kvec) == length(f̃max_kvec) == Nk
    @assert length(δΓ_kvec) == Nk

    for k in 1:Nk
        f̃_kvec[k] = constraint(S̃\SaddleVector(-d_kvec[k],zeros)) # Maybe move this to system constructor?
    end

    f̃ .= constraint(S̃\SaddleVector(-w,ψb))
    k_sheddingedges, activef̃limits_kvec = _findsheddingedges(Nk, e_kvec, f̃, f̃min_kvec, f̃max_kvec)

    # TODO: add loop below to correct vortex strengths (after adding constraint function to ConstrainedSystems)
    # while true
    #     f̃ .= constraint(S̃\SaddleVector(-w .+ sum(δΓ_kvec.*d_kvec,ψb))
    #     k_sheddingedges, activef̃limits_kvec = _findsheddingedges(Nk, e_kvec, f̃, f̃min_kvec, f̃max_kvec)
    #     if k_sheddingedges
    #         break
    #     end
    #     δΓ_kvec .= _computepointvortexstrengths(Nk, k_sheddingedges::Vector{Integer}, P_kvec, f̃_kvec, f̃lim_kvec, f₀, w)
    # end

    δΓ_kvec .= _computepointvortexstrengths(Nk, k_sheddingedges::Vector{Integer}, P_kvec, f̃_kvec, activef̃limits_kvec, f̃, f₀, w)

    _w_buf .= w .+ sum(δΓ_kvec.*d_kvec)

    steadyrhs = PotentialFlowRHS{T,TU,TF}(_w_buf, ψb, activef̃limits_kvec[k_sheddingedges], nothing, nothing)
    steadysys = RegularizedPotentialFlowSystem(Steady,S̃, f₀, e_kvec[k_sheddingedges], nothing, nothing, P_kvec[k_sheddingedges], zeros, ones, nothing)
    ldiv!(sol,steadysys,steadyrhs)

    return sol
end

function (\)(sys::UnregularizedPotentialFlowSystem,rhs::PotentialFlowRHS{T,TU,TF}) where {T,TU,TF}
    sol = PotentialFlowSolution(TU(),TF())
    ldiv!(sol,sys,rhs)
    return sol
end

function (\)(sys::RegularizedPotentialFlowSystem{FlowType,Nb,Nk,T,TU,TF,TE},rhs::PotentialFlowRHS{T,TU,TF}) where {FlowType,Nb,Nk,T,TU,TF,TE}
    if FlowType == Steady
        sol = PotentialFlowSolution(TU(),TF(),zeros(T,Nb),nothing)
    else
        sol = PotentialFlowSolution(TU(),TF(),zeros(T,Nb),zeros(T,Nk))
    end
    ldiv!(sol,sys,rhs)
    return sol
end

function _computesparsekuttaoperator(e::AbstractVector)::SparseMatrixCSC
    return sparse(I - ones(length(e))*e')
end

function _findsheddingedges(Nk, e_kvec, f̃, f̃min_kvec, f̃max_kvec)

    k_sheddingedges = Vector{Integer}()
    activef̃limits_kvec = Vector{Union{Nothing,eltype(f̃)}}()

    for k in 1:Nk
        if e_kvec[k]'*f̃ < f̃min_kvec[k]
            push!(activef̃limits_kvec,f̃min_kvec[k])
            push!(k_sheddingedges,k)
        elseif e_kvec[k]'*f̃ > f̃max_kvec[k]
            push!(activef̃limits_kvec,f̃max_kvec[k])
            push!(k_sheddingedges,k)
        else
            push!(activef̃limits_kvec,nothing)
        end
    end

    return k_sheddingedges, activef̃limits_kvec

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
