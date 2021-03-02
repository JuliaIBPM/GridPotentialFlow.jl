import LinearAlgebra: I, \, ldiv!, norm
import SparseArrays: AbstractSparseMatrix, SparseMatrixCSC, sparse
import Statistics: mean

export UnregularizedPotentialFlowSystem, RegularizedPotentialFlowSystem, PotentialFlowSystem, setd_kvec!

abstract type PotentialFlowSystem end

struct UnregularizedPotentialFlowSystem{Nb,T,TU,TF,TFB} <: PotentialFlowSystem
    S::SaddleSystem
    _TU_zeros::TU
    _TFB_ones::TFB # Matrix with the i,j-th entry equal to one if i is an index of the PointData that belongs to the j-th body and zero otherwise.
    _f_buf::TF

    function UnregularizedPotentialFlowSystem(S::SaddleSystem{T,Ns,Nc,TU,TF},_TF_ones::TFB) where {T,Ns,Nc,TU,TF,TFB}
        _TU_zeros = TU()
        _TU_zeros .= 0
        _f_buf = TF()
        new{size(_TF_ones,2),T,TU,TF,TFB}(S, _TU_zeros, _TF_ones, _f_buf)
    end
end

struct RegularizedPotentialFlowSystem{Nb,Nk,T,TU,TF,TE} <: PotentialFlowSystem
    S̃::SaddleSystem
    f₀::TF
    e_kvec::Vector{TE}
    d_kvec::Vector{TU}
    f̃_kvec::Vector{TF}
    P_kvec::Vector{AbstractSparseMatrix}
    _TF_zeros::TF
    _TF_ones::TF
    _w_buf::TU

    function RegularizedPotentialFlowSystem(S̃::SaddleSystem{T,Ns,Nc,TU,TF}, f₀::TF, e_kvec::Vector{TE}, d_kvec::Vector{TU}, f̃_kvec::Vector{TF}) where {T,Ns,Nc,TU,TF,TE}
        _w_buf = TU()
        P_kvec = _computesparsekuttaoperator.(e_kvec)
        _TF_zeros = TF()
        _TF_zeros .= 0
        _TF_ones = TF()
        _TF_ones .= 1
        new{1,length(e_kvec),T,TU,TF,TE}(S̃, f₀, e_kvec, d_kvec, f̃_kvec, P_kvec, _TF_zeros, _TF_ones, _w_buf)
    end

    # TODO check cost benefit of calling this constructor instead of the previous one because this one is unsafe
    function RegularizedPotentialFlowSystem(S̃::SaddleSystem{T,Ns,Nc,TU,TF}, f₀::TF, e_kvec::Vector{TE}, d_kvec::Vector{TU}, f̃_kvec::Vector{TF}, P_kvec::Vector{AbstractSparseMatrix}, _TF_zeros::TF, _TF_ones::TF, _w_buf::TU) where {T,Ns,Nc,TU,TF,TE}
        new{1,length(e_kvec),T,TU,TF,TE}(S̃, f₀, e_kvec, d_kvec, f̃_kvec, P_kvec, _TF_zeros, _TF_ones, _w_buf)
    end

end

function PotentialFlowSystem(S::SaddleSystem{T,Ns,Nc,TU,TF},_TF_ones) where {T,Ns,Nc,TU,TF}
    return UnregularizedPotentialFlowSystem(S,_TF_ones)
end

function PotentialFlowSystem(S̃::SaddleSystem{T,Ns,Nc,TU,TF}, f₀::TF, e_kvec::Vector{TE}) where {T,Ns,Nc,TU,TF,TE}
    return RegularizedPotentialFlowSystem(S̃,f₀,e_kvec,TU[],TF[])
end

function PotentialFlowSystem(S̃::SaddleSystem{T,Ns,Nc,TU,TF}, f₀::TF, e_kvec::Vector{TE}, d_kvec::Vector{TU}) where {T,Ns,Nc,TU,TF,TE}
    f̃_kvec = [TF() for i=1:length(e_kvec)]
    return RegularizedPotentialFlowSystem(S̃,f₀,e_kvec,d_kvec,f̃_kvec)
end

function ldiv!(sol::UnregularizedPotentialFlowSolution{T,TU,TF}, sys::UnregularizedPotentialFlowSystem{0,T,TU,TF}, rhs::UnregularizedPotentialFlowRHS{TU,TF}) where {T,TU,TF,TE}

    ldiv!(SaddleVector(sol.ψ,sol.f),sys.S,SaddleVector(-rhs.w,rhs.ψb))

    return sol
end

function ldiv!(sol::UnregularizedPotentialFlowSolution{T,TU,TF}, sys::UnregularizedPotentialFlowSystem{Nb,T,TU,TF}, rhs::UnregularizedPotentialFlowRHS{TU,TF}) where {Nb,T,TU,TF,TE}

    @unpack S, _TFB_ones, _TU_zeros, _f_buf = sys
    @unpack ψ, f, ψ₀ = sol
    @unpack w, ψb, Γb = rhs

    _computeconstraintonly!(f,S,ConstrainedSystems._unwrap_vec(-w),ψb,ConstrainedSystems._unwrap_vec(ψ))

    S₀ = Matrix(_TFB_ones'*(S.S⁻¹*_TFB_ones))
    ψ₀ = -S₀\(Γb-_TFB_ones'*f)

    ldiv!(SaddleVector(ψ,_f_buf),S,SaddleVector(_TU_zeros,_TFB_ones*ψ₀))

    f .= f .- _f_buf
    ψ .= reshape(-S.A⁻¹*(reshape(w,:) + S.B₁ᵀ*f),size(w))

    return sol
end

function ldiv!(sol::SteadyRegularizedPotentialFlowSolution{T,TU,TF}, sys::RegularizedPotentialFlowSystem{Nb,Nk,T,TU,TF,TE}, rhs::SteadyRegularizedPotentialFlowRHS{TU,TF,TSP}) where {Nb,Nk,T,TU,TF,TE,TSP}

    @unpack S̃, f₀, e_kvec, P_kvec, _TF_zeros, _TF_ones = sys
    @unpack ψ, f, ψ₀ = sol
    @unpack w, ψb, f̃lim_kvec = rhs

    # TODO: IMPLEMENT MULTIPLE BODIES

    @assert size(f̃lim_kvec) == size(e_kvec)

    _computeconstraintonly!(f,S̃,ConstrainedSystems._unwrap_vec(-w),ψb,ConstrainedSystems._unwrap_vec(ψ))
    ψ₀ .= mean(e_kvec)'*f .- mean(f̃lim_kvec)
    f .= mean(P_kvec)*f + _TF_ones*mean(f̃lim_kvec) # Represents f̃
    ψ .= reshape(-S̃.A⁻¹*(reshape(w,:) + S̃.B₁ᵀ*f),size(w))
    f .*= f₀

    return sol
end

function ldiv!(sol::UnsteadyRegularizedPotentialFlowSolution{T,TU,TF}, sys::RegularizedPotentialFlowSystem{Nb,Nk,T,TU,TF,TE}, rhs::UnsteadyRegularizedPotentialFlowRHS{TU,TF,TSP}) where {Nb,Nk,T,TU,TF,TE,TSP}

    @unpack S̃, f₀, e_kvec, d_kvec, f̃_kvec, P_kvec, _TF_zeros, _TF_ones, _w_buf = sys
    @unpack ψ, f, ψ₀, δΓ_kvec = sol
    @unpack w, ψb, f̃lim_kvec, Γw = rhs

    # TODO: IMPLEMENT MULTIPLE BODIES
    @assert length(d_kvec) == Nk
    @assert length(f̃_kvec) == Nk
    @assert length(f̃lim_kvec) == Nk
    @assert length(δΓ_kvec) == Nk

    _computeconstraintonly!(f,S̃,ConstrainedSystems._unwrap_vec(-w),ψb,ConstrainedSystems._unwrap_vec(ψ))

    activef̃lim_kvec = Vector{SuctionParameter}()
    for k in 1:Nk
        _computeconstraintonly!(f̃_kvec[k],S̃,ConstrainedSystems._unwrap_vec(-d_kvec[k]),_TF_zeros,ConstrainedSystems._unwrap_vec(ψ))
        activef̃lim = _findactivef̃limit(e_kvec[k],f,f̃lim_kvec[k])
        push!(activef̃lim_kvec,activef̃lim)
    end

    # The shedding edges are the _TF_ones for which the active f̃ limit is not Inf
    k_sheddingedges = [k for k in 1:Nk if activef̃lim_kvec[k] != Inf]

    # TODO: add loop below to correct vortex strengths (after adding constraint function to ConstrainedSystems)
    # while true
    #     f̃ .= constraint(S̃\SaddleVector(-w .+ sum(δΓ_kvec.*d_kvec,ψb))
    #     k_sheddingedges, activef̃limits_kvec = _findsheddingedges(Nk, e_kvec, f̃, f̃min_kvec, f̃max_kvec)
    #     if k_sheddingedges
    #         break
    #     end
    #     δΓ_kvec .= _computevortexstrengths(Nk, k_sheddingedges::Vector{Integer}, P_kvec, f̃_kvec, f̃lim_kvec, f₀, w)
    # end

    δΓ_kvec .= _computevortexstrengths(Nk, k_sheddingedges, P_kvec, f̃_kvec, activef̃lim_kvec, f, f₀, Γw)

    # Add the vorticity of the shedded vortices to the vorticity field and use this to compute the steady regularized solution
    _w_buf .= w .+ sum(δΓ_kvec.*d_kvec)

    # TODO: fails if there are no shedding edges
    steadyrhs = PotentialFlowRHS(_w_buf, ψb, activef̃lim_kvec[k_sheddingedges])
    steadysys = RegularizedPotentialFlowSystem(S̃, f₀, e_kvec[k_sheddingedges], TU[], TF[], P_kvec[k_sheddingedges], _TF_zeros, _TF_ones, _w_buf)
    steadysol = SteadyRegularizedPotentialFlowSolution(ψ,f,ψ₀)
    ldiv!(steadysol,steadysys,steadyrhs)

    return sol
end

function (\)(sys::UnregularizedPotentialFlowSystem, rhs::UnregularizedPotentialFlowRHS{TU,TF}) where {TU,TF}
    sol = PotentialFlowSolution(TU(),TF())
    ldiv!(sol,sys,rhs)
    return sol
end

function (\)(sys::RegularizedPotentialFlowSystem{Nb,Nk,T,TU,TF,TE}, rhs::SteadyRegularizedPotentialFlowRHS{TU,TF,TSP}) where {Nb,Nk,T,TU,TF,TE,TSP}
    sol = SteadyRegularizedPotentialFlowSolution(TU(),TF(),zeros(T,Nb))
    ldiv!(sol,sys,rhs)
    return sol
end

function (\)(sys::RegularizedPotentialFlowSystem{Nb,Nk,T,TU,TF,TE}, rhs::UnsteadyRegularizedPotentialFlowRHS{TU,TF,TSP}) where {Nb,Nk,T,TU,TF,TE,TSP}
    sol = UnsteadyRegularizedPotentialFlowSolution(TU(),TF(),zeros(T,Nb),zeros(T,Nk))
    ldiv!(sol,sys,rhs)
    return sol
end

function _computesparsekuttaoperator(e::AbstractVector)::SparseMatrixCSC
    return sparse(I - ones(length(e))*e')
end

function _findactivef̃limit(e::BodyUnitVector, f̃::TF, f̃lim::SuctionParameter) where {TF}

    f̃lim_range = SuctionParameterRange(-f̃lim,f̃lim)

    return _findactivef̃limit(e,f̃,f̃lim_range)
end

function _findactivef̃limit(e::BodyUnitVector, f̃::TF, f̃lim::SuctionParameterRange) where {TF}

    if e'*f̃ < f̃lim.min
        activef̃lim = f̃lim.min
    elseif e'*f̃ > f̃lim.max
        activef̃lim = f̃lim.max
    else
        activef̃lim = Inf
    end

    return activef̃lim
end

function _computevortexstrengths(Nk, k_sheddingedges::Vector{<:Integer}, P_kvec, f̃_kvec, f̃lim_kvec, f̃, f₀, Γw)

    δΓsys = zeros(Float64,Nk,Nk)
    δΓrhs = zeros(Float64,Nk)
    Γ₀ = sum(f₀)

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

# TODO: assign sys.d_kvec = d_kvec or remove d_kvec in vortexmodel altogether?
function setd_kvec!(sys::RegularizedPotentialFlowSystem{Nb,Nk,T,TU,TF,TE}, d_kvec::Vector{TU}) where {Nb,Nk,T,TU,TF,TE}

    resize!(sys.d_kvec,length(d_kvec))
    resize!(sys.f̃_kvec,length(d_kvec))
    sys.d_kvec .= d_kvec

end

function _removecirculation(sys,rhs)
    _TF_ones = typeof(rhs.ψb)()
    _TF_ones .= 1
    Γ₀ = _TF_ones'*(sys.S.S⁻¹*_TF_ones)
    rhs.ψb .= rhs.ψb-(_TF_ones*_TF_ones')/Γ₀*(sys.S.S⁻¹*rhs.ψb)
end

function _computeconstraintonly!(f,sys,w,ψb,nodes)
    nodes .= sys.A⁻¹*w

    sys.B₂A⁻¹r₁ .= sys.B₂*nodes
    sys._f_buf .= ψb
    sys._f_buf .-= sys.B₂A⁻¹r₁

    f .= sys.S⁻¹*sys._f_buf
end
