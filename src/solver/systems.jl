import LinearAlgebra: I, \, ldiv!, mul!, norm, LU, factorize, Factorization, dot
import SparseArrays: AbstractSparseMatrix, SparseMatrixCSC, sparse
import Statistics: mean

using LinearMaps

export UnregularizedPotentialFlowSystem, RegularizedPotentialFlowSystem, PotentialFlowSystem, setd_kvec!

abstract type AbstractPotentialFlowSystem{TU} end

abstract type PotentialFlowSystem end

struct IBPoisson{TU,TF} <: AbstractPotentialFlowSystem{TU}
    L::CartesianGrids.Laplacian
    R::RegularizationMatrix{TU,TF}
    E::InterpolationMatrix{TU,TF}
    _A⁻¹r₁::TU
    _B₁ᵀf::TU
    _f_buf::TF
    _Sfact::LU # Factorization

    function IBPoisson(L::CartesianGrids.Laplacian, R::RegularizationMatrix{TU,TF}, E::InterpolationMatrix{TU,TF}) where {TU,TF}
        _A⁻¹r₁ = TU()
        _B₁ᵀf = TU()
        _f_buf = TF()

        Rmap = LinearMap(x->vec(R*TF(x)),length(_A⁻¹r₁),length(_f_buf));
        Emap = LinearMap(x->vec(E*TU(x)),length(_f_buf),length(_A⁻¹r₁))
        L⁻¹map = LinearMap(x->vec(L\TU(x)),length(_A⁻¹r₁));
        Smap = -Emap*L⁻¹map*Rmap
        _Sfact = factorize(Matrix(Smap))

        new{TU,TF}(L, R, E, _A⁻¹r₁, _B₁ᵀf, _f_buf, _Sfact)
    end
end

"""
$(TYPEDEF)

... System with a single type of constraints on f
"""
struct ConstrainedIBPoisson{TU,TF,TB1T_2,TB2_2} <: AbstractPotentialFlowSystem{TU}
    ibp::IBPoisson{TU,TF}
    B₁ᵀpart2::TB1T_2
    B₂part2::TB2_2
    _f_buf::TF
    _f_buf2::TF
    _S₀fact::Union{Factorization,Float64}
    _ψ₀_buf::Vector{Float64}
    _sol_buf::IBPoissonSolution{TU,TF}

    function ConstrainedIBPoisson(L::CartesianGrids.Laplacian, R::RegularizationMatrix{TU,TF}, E::InterpolationMatrix{TU,TF}, B₁ᵀpart2::TB1T_2, B₂part2::TB2_2) where {TU,TF,TB1T_2,TB2_2}
        Nb = size(B₁ᵀpart2,2) # number of bodies
        ibp = IBPoisson(L,R,E)
        S₀ = Matrix{Float64}(undef,Nb,Nb)
        _S₀fact = _computeS₀fact(S₀,ibp._Sfact,B₁ᵀpart2,B₂part2,deepcopy(B₁ᵀpart2))
        _ψ₀_buf = zeros(Nb)
        _sol_buf = IBPoissonSolution{TU,TF}(TU(),TF())
        new{TU,TF,TB1T_2,TB2_2}(ibp, B₁ᵀpart2, B₂part2, TF(), TF(), _S₀fact, _ψ₀_buf, _sol_buf)
    end
end

struct SteadyRegularizedIBPoisson{TU,TF,TB1T_2,TB2_2} <: AbstractPotentialFlowSystem{TU}
    cibp::ConstrainedIBPoisson{TU,TF,TB1T_2,TB2_2}
    f₀::TF

    function SteadyRegularizedIBPoisson(L::CartesianGrids.Laplacian, R::RegularizationMatrix{TU,TF}, E::InterpolationMatrix{TU,TF}, B₁ᵀpart2::TB1T_2, B₂part2::TB2_2) where {TU,TF,TB1T_2,TB2_2}

        _TU_buf = TU()

        f₀ = TF()
        _TF_buf = TF()
        _TF_buf .= 1.0
        ibp = IBPoisson(L, R, E)
        ldiv!(IBPoissonSolution(_TU_buf,f₀),ibp,IBPoissonRHS(_TU_buf,_TF_buf), onlyf=true, zerow=true)

        R̃ = deepcopy(R)
        R̃.M .= R̃.M*Diagonal(f₀)
        cibp = ConstrainedIBPoisson(L, R̃, E, B₁ᵀpart2, B₂part2)
        new{TU,TF,TB1T_2,TB2_2}(cibp, f₀)
    end
end

struct UnsteadyRegularizedIBPoisson{Nb,Ne,TU,TF} <: AbstractPotentialFlowSystem{TU}
    ibp::IBPoisson{TU,TF}
    one_vec::Vector{TF}
    e_vec::Vector{TF}
    f₀::TF
    f₀_vec::Vector{TF}
    f̃₀_vec::Vector{TF}
    f̃_vec::Vector{TF}
    Γ₀::Float64
    d_vec::Vector{TU}
    _activef̃lim_vec::Vector{f̃Limit}
    _TU_buf::TU
    _f_buf::TF
    _Souter::Matrix{Float64}
    _r₂_buf::Vector{Float64}
    _y_buf::Vector{Float64}

    function UnsteadyRegularizedIBPoisson(L::CartesianGrids.Laplacian, R::RegularizationMatrix{TU,TF}, E::InterpolationMatrix{TU,TF}, one_vec::Vector{TF}, e_vec::Vector{TF}) where {TU,TF}
        Nb = length(one_vec) # number of bodies
        Ne = length(e_vec) # number of edges

        ibp = IBPoisson(L, R, E)

        _TU_buf = TU()
        _f_buf = TF()
        _f_buf2 = TF()

        f₀ = TF()

        _f_buf .= 1.0
        ldiv!(IBPoissonSolution(_TU_buf,f₀),ibp,IBPoissonRHS(_TU_buf,_f_buf), onlyf=true, zerow=true)

        f₀_vec = [TF() for i=1:Nb]
        for b in 1:Nb
            ldiv!(IBPoissonSolution(_TU_buf,f₀_vec[b]),ibp,IBPoissonRHS(_TU_buf,one_vec[b]), onlyf=true, zerow=true)
        end

        R̃ = deepcopy(R)
        R̃.M .= R̃.M*Diagonal(f₀)
        ibp = IBPoisson(L, R̃, E)

        f̃₀_vec = [TF() for i=1:Nb]
        for b in 1:Nb
            ldiv!(IBPoissonSolution(_TU_buf,f̃₀_vec[b]),ibp,IBPoissonRHS(_TU_buf,one_vec[b]), onlyf=true, zerow=true)
        end

        d_vec = [TU() for i=1:Ne]
        f̃_vec = [TF() for i=1:Ne]
        Γ₀ = sum(f₀)
        _activef̃lim_vec = zeros(Ne)
        _Souter = zeros(Nb+Ne,Nb+Ne)
        _r₂_buf = zeros(Nb+Ne)
        _y_buf = zeros(Nb+Ne)

        new{Nb,Ne,TU,TF}(ibp, one_vec, e_vec, f₀, f₀_vec, f̃₀_vec, f̃_vec, Γ₀, d_vec, _activef̃lim_vec, _TU_buf, _f_buf, _Souter, _r₂_buf, _y_buf)
    end
end

function ldiv!(sol::TS, sys::IBPoisson{TU,TF}, rhs::TR; zerow=false, zeroψb=false, useprovidedf=false, onlyf=false) where {TU,TF,TS,TR}

    if useprovidedf # in this case we use the f provided in sol.f and only calculate ψ = L⁻¹(-w-Rf)
        mul!(sys._B₁ᵀf, sys.R, sol.f)
        sys._B₁ᵀf .= .-rhs.w .- sys._B₁ᵀf
        ldiv!(sol.ψ, sys.L, sys._B₁ᵀf)
    else
        if zerow
            sys._f_buf .= 0.0
        else
            ldiv!(sys._A⁻¹r₁, sys.L, rhs.w) # -A⁻¹r₁ (note minus sign because w and not -w is used)
            mul!(sys._f_buf, sys.E, sys._A⁻¹r₁) # -B₂A⁻¹r₁ (note minus sign)
        end

        if !zeroψb
            sys._f_buf .= rhs.ψb .+ sys._f_buf # r₂ - B₂A⁻¹r₁
        end

        ldiv!(sol.f, sys._Sfact, sys._f_buf) # S⁻¹(r₂ - B₂A⁻¹r₁)

        if !onlyf
            mul!(sys._B₁ᵀf, sys.R, sol.f)
            ldiv!(sol.ψ, sys.L, sys._B₁ᵀf)
            sol.ψ .= sys._A⁻¹r₁ .- sol.ψ
        end
    end

    return sol
end

function ldiv!(sol::TS, sys::ConstrainedIBPoisson{TU,TF,TB1T_2,TB2_2}, rhs::ConstrainedIBPoissonRHS) where {TU,TF,TB1T_2,TB2_2,TS}

    # fstar = S⁻¹(ψb+EL⁻¹)
    ldiv!(sol, sys.ibp, rhs, onlyf=true)

    # ψ₀ = S₀\(Γb-B₂₂*fstar)
    mul!(sys._ψ₀_buf, sys.B₂part2, sol.f)
    sys._ψ₀_buf .= rhs.fconstraintRHS .- sys._ψ₀_buf
    ldiv!(sol.ψ₀, sys._S₀fact, sys._ψ₀_buf)

    # f = fstar - S⁻¹(B₁ᵀ₂*ψ₀)
    mul!(sys._f_buf, sys.B₁ᵀpart2, sol.ψ₀)
    ldiv!(sys._f_buf2, sys.ibp._Sfact, sys._f_buf)
    sol.f .-= sys._f_buf2

    # ψ = L⁻¹(-w - Rf)
    ldiv!(sol, sys.ibp, rhs, useprovidedf=true)

end

function ldiv!(sol::TS, sys::SteadyRegularizedIBPoisson{TU,TF,TB1T_2,TB2_2}, rhs::TR) where {TU,TF,TB1T_2,TB2_2,TS,TR}

    ldiv!(sol, sys.cibp, rhs)
    sol.f .*= sys.f₀

end

function ldiv!(sol::TS, sys::UnsteadyRegularizedIBPoisson{Nb,Ne,TU,TF}, rhs::TR) where {Nb,Ne,TU,TF,TS,TR}

    ldiv!(sol, sys.ibp, rhs, onlyf=true)

    for k in 1:Ne
        # eq 2.47
        ldiv!(IBPoissonSolution(sys._TU_buf, sys.f̃_vec[k]), sys.ibp, IBPoissonRHS(sys.d_vec[k],sys._f_buf), zeroψb=true, onlyf=true)
        # eq 2.60
        sys._activef̃lim_vec[k] = _findactivef̃limit(sys.e_vec[k], sol.f, rhs.f̃lim_vec[k])
    end

    # The shedding edges are the ones for which the active f̃ limit is not Inf
    sheddingedges = [e for e in 1:length(sys._activef̃lim_vec) if sys._activef̃lim_vec[e] != Inf]

    # eq 2.49, eq 2.56, eq 2.61
    _computeδΓandψ₀!(sol,sys,rhs,sheddingedges)

    # f̃ = f̃star + Σf̃ᵢδΓᵢ - Σf̃₀ⱼψ₀ⱼ
    _addlincomboto!(sol.f, sol.δΓ_vec, sys.f̃_vec)
    sys._f_buf .= 0.0
    _addlincomboto!(sys._f_buf, sol.ψ₀, sys.f̃₀_vec)
    sol.f .-= sys._f_buf

    # ψ = L⁻¹(-w - ΣdᵢδΓᵢ - Rf)
    sys._TU_buf .= rhs.w
    _addlincomboto!(sys._TU_buf, sol.δΓ_vec, sys.d_vec)
    ldiv!(sol, sys.ibp, IBPoissonRHS(sys._TU_buf, sys._f_buf), useprovidedf=true)

    sol.f .*= sys.f₀

end

function _addlincomboto!(s,c,A)
    for i in 1:length(c)
        s .+= c[i].*A[i]
    end
end

# function _computeS₀fact(S₀,Sfact,B₁ᵀ₂cols,B₂₁rows,_S⁻¹B₁ᵀ₂)
#     for c in 1:length(B₁ᵀ₂cols)
#         ldiv!(view(_S⁻¹B₁ᵀ₂,:,c),Sfact,B₁ᵀ₂cols[c])
#     end
#
#     for r in 1:length(B₂₁rows)
#         mul!(view(S₀,r,:),B₂₁rows[r],TB1T_2_buf)
#     end
#
#     S₀fact = factorize(S₀)
# end

function _computeS₀fact(S₀,Sfact,B₁ᵀpart2,B₂part2,TB1T_2_buf)
    ldiv!(TB1T_2_buf,Sfact,B₁ᵀpart2)
    mul!(S₀,B₂part2,TB1T_2_buf)
    S₀ .= .-S₀
    S₀fact = factorize(S₀)
end

function ldiv!(x::ScalarData ,A::LU, f::ScalarData)
    ldiv!(x.data,A,f.data)
end

function _computeδΓandψ₀!(sol::TS, sys::UnsteadyRegularizedIBPoisson{Nb,Ne,TU,TF}, rhs::TR, sheddingedges::Vector{Int}) where {Nb,Ne,TU,TF,TS,TR}
    sys._Souter .= 0.0
    for i in 1:Ne
        if i in sheddingedges
            for b in 1:Nb
                sys._Souter[i,b] = dot(sys.e_vec[i], sys.f̃₀_vec[b])
            end
            for e in 1:Ne
                sys._Souter[i,end-Ne+e] = -dot(sys.e_vec[i], sys.f̃_vec[e])
            end
            sys._r₂_buf[i] = sys._activef̃lim_vec[i] - dot(sys.e_vec[i], sol.f)
        else
            sys._Souter[end-Ne+i] = 1.0
            sys._r₂_buf[i] = 0.0
        end
    end
    for i in 1:Nb
        for b in 1:Nb
            sys._Souter[end-Nb+i,b] = dot(sys.f₀_vec[i], sys.f̃₀_vec[b])
        end
        for e in 1:Ne
            sys._Souter[end-Nb+i,end-Ne+e] = -dot(sys.f₀_vec[i], sys.f̃_vec[e])
        end
        sys._r₂_buf[end-Nb+i] = -rhs.Γw_vec[i] - dot(sys.f₀_vec[i], sol.f)
    end
    sys._Souter .= .-sys._Souter
    for i in 1:Nb
        sys._Souter[end-Nb+i,end-Nb+i] += 1.0
    end

    ldiv!(sys._y_buf, factorize(sys._Souter), sys._r₂_buf)
    sol.ψ₀ .= sys._y_buf[1:Ne]
    sol.δΓ_vec .= sys._y_buf[end-Nb+1:end]
end













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
        # _TFTFt_ones = _TF_ones*_TF_ones'
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

function ldiv!(sol::SteadyRegularizedPotentialFlowSolution{T,TU,TF}, sys::RegularizedPotentialFlowSystem{Nb,Nk,T,TU,TF,TE}, rhs::SteadyRegularizedPotentialFlowRHS{TU,TF}) where {Nb,Nk,T,TU,TF,TE}

    @unpack S̃, f₀, e_kvec, P_kvec, _TF_zeros, _TF_ones, _w_buf = sys
    @unpack ψ, f, ψ₀ = sol
    @unpack w, ψb, f̃lim_kvec = rhs

    # TODO: IMPLEMENT MULTIPLE BODIES

    @assert size(f̃lim_kvec) == size(e_kvec)

    # eq 2.33
    _computeconstraintonly!(f,S̃,ConstrainedSystems._unwrap_vec(-w),ψb,ConstrainedSystems._unwrap_vec(ψ))
    # eq 2.35 + accounting for generalized edge conditions
    ψ₀ .= mean(e_kvec)'*f .- mean(f̃lim_kvec)
    # eq 2.37 + accounting for generalized edge conditions
    f .= mean(P_kvec)*f + _TF_ones*mean(f̃lim_kvec) # Represents f̃
    # eq 2.38
    ψ .= reshape(reshape(_w_buf,:) .+ S̃.B₁ᵀ*f,size(ψ))
    ψ .= reshape(-S̃.A⁻¹*reshape(ψ,:),size(ψ))
    f .*= f₀

    return sol
end

function ldiv!(sol::UnsteadyRegularizedPotentialFlowSolution{T,TU,TF}, sys::RegularizedPotentialFlowSystem{Nb,Nk,T,TU,TF,TE}, rhs::UnsteadyRegularizedPotentialFlowRHS{TU,TF}) where {Nb,Nk,T,TU,TF,TE}

    @unpack S̃, f₀, e_kvec, d_kvec, f̃_kvec, P_kvec, _TF_zeros, _TF_ones, _w_buf = sys
    @unpack ψ, f, ψ₀, δΓ_kvec = sol
    @unpack w, ψb, f̃lim_kvec, Γw = rhs

    # TODO: IMPLEMENT MULTIPLE BODIES
    @assert length(d_kvec) == Nk
    @assert length(f̃_kvec) == Nk
    @assert length(f̃lim_kvec) == Nk
    @assert length(δΓ_kvec) == Nk

    # eq 2.33
    _computeconstraintonly!(f,S̃,ConstrainedSystems._unwrap_vec(-w),ψb,ConstrainedSystems._unwrap_vec(ψ))

    activef̃lim_kvec = Vector{SuctionParameter}()
    for k in 1:Nk
        # eq 2.47
        _computeconstraintonly!(f̃_kvec[k],S̃,ConstrainedSystems._unwrap_vec(-d_kvec[k]),_TF_zeros,ConstrainedSystems._unwrap_vec(ψ))
        activef̃lim = _findactivef̃limit(e_kvec[k],f,f̃lim_kvec[k])
        push!(activef̃lim_kvec,activef̃lim)
    end

    # The shedding edges are the _TF_ones for which the active f̃ limit is not Inf
    k_sheddingedges = [k for k in 1:Nk if activef̃lim_kvec[k] != Inf]

    # eq 2.49, eq 2.56, eq 2.61
    δΓ_kvec .= _computevortexstrengths(Nk, k_sheddingedges, P_kvec, f̃_kvec, activef̃lim_kvec, f, f₀, Γw)

    # Add the vorticity of the shedded vortices to the vorticity field and use this to compute the steady regularized solution
    _w_buf .= w .+ sum(δΓ_kvec.*d_kvec)
    # Add the bound vortex sheet strength induced by the shedded vortices to the bound vortex sheet strength induced by w and the boundary conditions
    f .= f .+ sum(δΓ_kvec.*f̃_kvec)
    # eq 2.57 + accounting for generalized edge conditions
    ψ₀ .= mean(e_kvec[k_sheddingedges])'*f .- mean(activef̃lim_kvec[k_sheddingedges])
    # eq 2.58 + accounting for generalized edge conditions
    f .= mean(P_kvec[k_sheddingedges])*f + _TF_ones*mean(activef̃lim_kvec[k_sheddingedges]) # Represents f̃
    ψ .= reshape(reshape(_w_buf,:) .+ S̃.B₁ᵀ*f,size(ψ))
    ψ .= reshape(-S̃.A⁻¹*reshape(ψ,:),size(ψ))
    # eq 2.24
    f .*= f₀

    return sol
end

function (\)(sys::UnregularizedPotentialFlowSystem, rhs::UnregularizedPotentialFlowRHS{TU,TF}) where {TU,TF}
    sol = PotentialFlowSolution(TU(),TF())
    ldiv!(sol,sys,rhs)
    return sol
end

function (\)(sys::RegularizedPotentialFlowSystem{Nb,Nk,T,TU,TF,TE}, rhs::SteadyRegularizedPotentialFlowRHS{TU,TF}) where {Nb,Nk,T,TU,TF,TE}
    sol = SteadyRegularizedPotentialFlowSolution(TU(),TF(),zeros(T,Nb))
    ldiv!(sol,sys,rhs)
    return sol
end

function (\)(sys::RegularizedPotentialFlowSystem{Nb,Nk,T,TU,TF,TE}, rhs::UnsteadyRegularizedPotentialFlowRHS{TU,TF}) where {Nb,Nk,T,TU,TF,TE}
    sol = UnsteadyRegularizedPotentialFlowSolution(TU(),TF(),zeros(T,Nb),zeros(T,Nk))
    ldiv!(sol,sys,rhs)
    return sol
end

function _computesparsekuttaoperator(e::AbstractVector)::SparseMatrixCSC
    return sparse(I - ones(length(e))*e')
end

function _findactivef̃limit(e::BodyUnitVector, f̃::TF, f̃lim::f̃Limits) where {TF}

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
function setd_kvec!(sys::RegularizedPotentialFlowSystem{Nb,Nk,T,TU,TF,TE}, d_kvec::Vector) where {Nb,Nk,T,TU,TF,TE}
    @assert size(d_kvec) == size(sys.d_kvec)

    resize!(sys.d_kvec,length(d_kvec))
    resize!(sys.f̃_kvec,length(d_kvec))
    sys.d_kvec .= d_kvec

end

# function _removecirculation!(ψb,sys)
#     Γ₀ = sum(sys.f₀)
#     ψb .= ψb .- 1/Γ₀*sys._TFTFt_ones*sys.S̃.S⁻¹*ψb
# end

function _computeconstraintonly!(f,sys,w,ψb,nodes)
    nodes .= sys.A⁻¹*w

    sys.B₂A⁻¹r₁ .= sys.B₂*nodes
    sys._f_buf .= ψb
    sys._f_buf .-= sys.B₂A⁻¹r₁

    f .= sys.S⁻¹*sys._f_buf
end






# function ldiv!(sol::TS, sys::ConstrainedIBPoisson{TU,TF,TB1T_2,TB2_2}, rhs::TR; useprovidedψf=false, useprovidedψ₀=false) where {TU,TF,TB1T_2,TB2_2,TS,TR}
#
#     # USE RHS GET METHOD
#     if rhs isa UnregularizedPotentialFlowRHS
#         fconstraintRHS = rhs.Γb
#     elseif rhs isa SteadyRegularizedPotentialFlowRHS
#         fconstraintRHS = rhs.f̃lim_kvec
#     end
#
#     # _computeconstraintonly!(f,S,ConstrainedSystems._unwrap_vec(-w),ψb,ConstrainedSystems._unwrap_vec(ψ))
#     if !useprovidedψf
#         ldiv!(sol, sys.ibp, rhs)
#     end
#
#     # ψ₀ = -S₀\(Γb-B₂*f)
#     if !useprovidedψ₀
#         mul!(sys._ψ₀_buf, sys.B₂part2, sol.f)
#         sys._ψ₀_buf .= fconstraintRHS .- sys._ψ₀_buf
#         ldiv!(sol.ψ₀, sys._S₀fact, sys._ψ₀_buf)
#         sol.ψ₀ .= .-sol.ψ₀
#     end
#
#     # ldiv!(SaddleVector(ψ,_f_buf),S,SaddleVector(_TU_zeros,_TB2_2_ones*ψ₀))
#     mul!(sys._f_buf, sys.B₁ᵀpart2, sol.ψ₀)
#     ldiv!(sys._sol_buf, sys.ibp, IBPoissonRHS(rhs.w, sys._f_buf), zerow=true)
#
#     sol.f .= sol.f .- sys._sol_buf.f
#     sol.ψ .= sol.ψ .- sys._sol_buf.ψ
#
# end
