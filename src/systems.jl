import LinearAlgebra: I, \, ldiv!
import SparseArrays: AbstractSparseMatrix, SparseMatrixCSC, sparse
import Statistics: mean

abstract type FlowType end
abstract type Steady <: FlowType end
abstract type Unsteady <: FlowType end

struct UnregularizedPotentialFlowSystem
    S::SaddleSystem
end

struct RegularizedPotentialFlowSystem{FlowType,T,TU,TF,TE}
    S̃::SaddleSystem
    f₀::TF
    e_kvec::Vector{TE}
    d_kvec::Union{Nothing,Vector{TU}}
    f̃_kvec::Union{Nothing,Vector{TF}}
    P_kvec::Vector{AbstractSparseMatrix}
    zeros::TF
    ones::TF

    function RegularizedPotentialFlowSystem(flowtype::Type{<:FlowType},S̃::SaddleSystem{T,Ns,Nc,TU,TF}, f₀::TF, e_kvec::Vector{TE}, d_kvec::Union{Nothing,Vector{TU}}=nothing, f̃_kvec::Union{Nothing,Vector{TF}}=nothing) where {T,Ns,Nc,TU,TF,TE}
        if flowtype == Unsteady
            if isnothing(d_kvec) error("d_kvec is not defined") end
            if isnothing(f̃_kvec) error("f̃_kvec is not defined") end
            if size(e_kvec) != size(d_kvec) error("Incompatible number of elements in e_kvec and d_kvec") end
            if size(e_kvec) != size(f̃_kvec) error("Incompatible number of elements in e_kvec and f̃_kvec") end
        end
        P_kvec = _computesparsekuttaoperator.(e_kvec)
        zeros = TF()
        zeros .= 0
        ones = TF()
        ones .= 1
        new{flowtype,T,TU,TF,TE}(S̃, f₀, e_kvec, d_kvec, f̃_kvec, P_kvec, zeros, ones)
    end

    # TODO check cost benefit of calling this constructor instead of the previous one because this one is unsafe
    function RegularizedPotentialFlowSystem(flowtype::Type{<:FlowType},S̃::SaddleSystem{T,Ns,Nc,TU,TF}, f₀::TF, e_kvec::Vector{TE}, d_kvec::Union{Nothing,Vector{TU}}, f̃_kvec::Union{Nothing,Vector{TF}}, P_kvec::Vector{AbstractSparseMatrix}, zeros::TF, ones::TF) where {T,Ns,Nc,TU,TF,TE}
        new{flowtype,T,TU,TF,TE}(S̃, f₀, e_kvec, d_kvec, f̃_kvec, P_kvec, zeros, ones)
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

function ldiv!(sol::PotentialFlowSolution{T,TU,TF,TB}, sys::RegularizedPotentialFlowSystem{Steady,T,TU,TF,TE}, rhs::PotentialFlowRHS{T,TU,TF}) where {FlowType,T,TU,TF,TE,TB}

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

function ldiv!(sol::PotentialFlowSolution{T,TU,TF,TB}, sys::RegularizedPotentialFlowSystem{Unsteady,T,TU,TF,TE}, rhs::PotentialFlowRHS{T,TU,TF}) where {FlowType,T,TU,TF,TE,TB}

    @unpack S̃, f₀, e_kvec, d_kvec, f̃_kvec, P_kvec, zeros, ones = sys
    @unpack ψ, f̃, ψ₀, δΓ_kvec = sol
    @unpack w, ψb, f̃min_kvec, f̃max_kvec = rhs

    # TODO: IMPLEMENT MULTIPLE BODIES

    @assert size(f̃min_kvec) == size(f̃max_kvec) == size(e_kvec)
    @assert size(δΓ_kvec) == size(e_kvec)

    for k in 1:length(e_kvec)
        f̃_kvec[k] .= constraint(S̃\SaddleVector(-d_kvec[k],zeros)) # Maybe move this to system constructor?
    end

    f̃ .= constraint(S̃\SaddleVector(-w,ψb))
    k_sheddingedges, activef̃limits_kvec = _findsheddingedges(e_kvec, f̃, f̃min_kvec, f̃max_kvec)

    # while true
    #     f̃ .= constraint(S̃\SaddleVector(-w .+ sum(δΓ_kvec.*d_kvec,ψb))
    #     k_sheddingedges, activef̃limits_kvec = _findsheddingedges(e_kvec, f̃, f̃min_kvec, f̃max_kvec)
    #     if k_sheddingedges
    #         break
    #     end
    #     δΓ_kvec .= _computepointvortexstrengths(k_sheddingedges::Vector{Integer}, e_kvec, P_kvec, f̃_kvec, f̃lim_kvec, f₀, w)
    # end

    println(k_sheddingedges)
    println(activef̃limits_kvec)

    δΓ_kvec .= _computepointvortexstrengths(k_sheddingedges::Vector{Integer}, e_kvec, P_kvec, f̃_kvec, activef̃limits_kvec, f̃, f₀, w)

    w .= w .+ sum(δΓ_kvec.*d_kvec)

    steadyrhs = PotentialFlowRHS{T,TU,TF}(w, ψb, activef̃limits_kvec[k_sheddingedges], nothing, nothing)
    steadysys = RegularizedPotentialFlowSystem(Steady,S̃, f₀, e_kvec[k_sheddingedges], nothing, nothing, P_kvec[k_sheddingedges], zeros, ones)
    ldiv!(sol,steadysys,steadyrhs)

    return sol
end

# function ldiv!(sol::PotentialFlowSolution{T,TU,TF}, sys::RegularizedPotentialFlowSystem{FlowType,T,TU,TF,TE}, rhs::PotentialFlowRHS{T,TU,TF}) where {FlowType,T,TU,TF,TE}
#
#     @unpack S̃, f₀, e_kvec, d_kvec, f̃_kvec, P_kvec, zeros, ones = sys
#     @unpack ψ, f̃, ψ₀ = sol
#     @unpack w, ψb = rhs
#
#     # TODO: maybe type check here instead of in arguments list?
#     # TODO: IMPLEMENT MULTIPLE BODIES
#
#
#     # Streamfunction and vortex sheet strength without shedded vortices
#     ψ,f̃ = S̃\(-w,ψb);
#
#     if FlowType == Unsteady
#         δΓ_kvec, f̃_kvec, f̃activelimit_kvec = _computepointvortexstrengths!(sol,sys,rhs)
#         f̃ = f̃ + sum(δΓ_kvec.*f̃_kvec) # Maybe remove this one and put next line in before ψ,f̃ = S̃\(-w,ψb)
#         w = w + sum(δΓ_kvec.*d_kvec)
#     else
#         @assert size(rhs.f̃limit_kvec
#         f̃activelimit_kvec = rhs.f̃limit_kvec
#     end
#
#     ψ₀ = mean(e_kvec)'*f̃ - mean(f̃activelimit_kvec); # is ψ₀
#     f̃ .= mean(P_kvec)*f̃ + ones*mean(f̃activelimit_kvec);
#     ψ .= reshape(-S̃.A⁻¹*(reshape(w,:) + S̃.B₁ᵀ*f̃),size(ψ))
#
#     return sol
# end

# function (\)(sys::PotentialFlowSystem{Nk,T},rhs::Tuple{AbstractMatrix, AbstractVector, AbstractVector, Real}) where {Nk,T}
#     negw, ψb, f̃limk, negΓw = rhs
#     # ONLY WORKS FOR A SINGLE BODY FOR NOW!
#     sol = (similar(negw),similar(ψb),Array{T}(undef, 1),Array{T}(undef, Nk))
#     ldiv!(sol,sys,rhs)
#     return sol
# end

function _computesparsekuttaoperator(e::AbstractVector)::SparseMatrixCSC
    return sparse(I - ones(length(e))*e')
end

function _findsheddingedges(e_kvec, f̃, f̃min_kvec, f̃max_kvec)

    k_sheddingedges = Vector{Integer}()
    activef̃limits_kvec = Vector{Union{Nothing,eltype(f̃)}}()

    for k in 1:length(e_kvec)
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

function _computepointvortexstrengths(k_sheddingedges::Vector{Integer}, e_kvec, P_kvec, f̃_kvec, f̃lim_kvec, f̃, f₀, w)
    Nk = length(e_kvec)

    δΓsys = zeros(eltype(w),Nk,Nk)
    δΓrhs = zeros(eltype(w),Nk)
    Γ₀ = sum(f₀)
    Γw = sum(w)

    for k in 1:length(e_kvec)
        if k in k_sheddingedges
            δΓsys[k_sheddingedges,k] .= 1 + f₀'*P_kvec[k]*f̃_kvec[k]
            δΓrhs[k] = Γw + f₀'*P_kvec[k]*f̃ + Γ₀*f̃lim_kvec[k]
        else
            δΓsys[k,k] = 1
            δΓrhs[k] = 0
        end
    end

    δΓ_kvec = -δΓsys\δΓrhs

    println(δΓsys)
    println(δΓrhs)
    println(δΓ_kvec)

    return δΓ_kvec
end

# function _computepointvortexstrengths!(sol::PotentialFlowSolution{T,TU,TF}, sys::RegularizedPotentialFlowSystem{FlowType,T,TU,TF,TE}, rhs::PotentialFlowRHS{T,TU,TF}) where {T,TU,TF,TE}
#     @unpack S̃, f₀, e_kvec, d_kvec, f̃_kvec, P_kvec, zeros, ones = sys
#     @unpack ψ, f̃, ψ₀, δΓ_kvec = sol
#     @unpack w, ψb, f̃min_kvec, f̃max_kvec, f̃lim_kvec = rhs
#
#     f̃lim_kvec = fill(0.0,length(e_kvec))
#     releasevortex = fill(Bool(0),length(e_kvec))
#
#     δΓsysrow = T[]
#     δΓrhs = T[]
#     indicestocompute = Integer[]
#     Γ₀ = sum(f₀)
#     Γw = sum(w)
#
#     for i in 1:length(e_kvec)
#         f̃_kvec[i] .= constraint(S̃\SaddleVector(-d_kvec[i],zeros));
#         f̃min_kvec[i] ≤ e_kvec[i]'*f̃ ≤ f̃max_kvec[i] ? (releasevortex[i] = false; δΓ_kvec[i] = 0.0) : (releasevortex[i] = true; e_kvec[i]'*f̃ < f̃min_kvec[i] ? f̃lim_kvec[i] = f̃min_kvec[i] : f̃lim_kvec[i] = f̃max_kvec[i])
#         if releasevortex[i]
#             push!(indicestocompute,i)
#             push!(δΓsysrow, 1 + f₀'*P_kvec[i]*f̃_kvec[i])
#             push!(δΓrhs, Γw + f₀'*P_kvec[i]*f̃ + Γ₀*f̃lim_kvec[i])
#         end
#     end
#
#     δΓsys = repeat(δΓsysrow',length(indicestocompute))
#     δΓ_kvec[indicestocompute] .= -δΓsys\δΓrhs
#
#     return δΓ_kvec, f̃_kvec, f̃lim_kvec
# end
#
# function _computepointvortexstrengths!(sol::PotentialFlowSolution{T,TU,TF}, sys::RegularizedPotentialFlowSystem{FlowType,T,TU,TF,TE}, rhs::PotentialFlowRHS{T,TU,TF}) where {T,TU,TF,TE}
#     @unpack S̃, f₀, e_kvec, d_kvec, f̃_kvec, P_kvec, zeros, ones = sys
#     @unpack ψ, f̃, ψ₀, δΓ_kvec = sol
#     @unpack w, ψb, f̃min_kvec, f̃max_kvec, f̃lim_kvec = rhs
#
#     Nk = length(e_kvec)
#
#     δΓsys = zeros(T,Nk,Nk)
#     δΓrhs = zeros(T,Nk)
#     Γ₀ = sum(f₀)
#     Γw = sum(w)
#
#     for i in 1:length(e_kvec)
#         if f̃min_kvec[i] ≤ e_kvec[i]'*f̃ ≤ f̃max_kvec[i]
#             δΓsys[i,i] = 0
#             δΓrhs[i] = 0
#         else
#             if e_kvec[i]'*f̃ < f̃min_kvec[i]
#                 f̃lim_kvec[i] = f̃min_kvec[i]
#             else
#                 f̃lim_kvec[i] = f̃max_kvec[i]
#             end
#             δΓsys[:,i] .= 1 + f₀'*P_kvec[i]*f̃_kvec[i]
#             δΓrhs[i] = Γw + f₀'*P_kvec[i]*f̃ + Γ₀*f̃lim_kvec[i]
#         end
#     end
#
#     δΓ_kvec .= -δΓsys\δΓrhs
#
#     return δΓ_kvec, f̃lim_kvec
# end
