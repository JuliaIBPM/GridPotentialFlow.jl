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
struct ConstrainedIBPoisson{Nb,TU,TF} <: AbstractPotentialFlowSystem{TU}
    ibp::IBPoisson{TU,TF}
    B₁ᵀ₂cols::Vector{TF}
    B₂₂rows::Vector{TF}
    _f_buf::TF
    _f_buf2::TF
    _S₀fact::Union{Factorization,Float64}
    _ψ₀_buf::Vector{Float64}
    _sol_buf::IBPoissonSolution{TU,TF}

    function ConstrainedIBPoisson(L::CartesianGrids.Laplacian, R::RegularizationMatrix{TU,TF}, E::InterpolationMatrix{TU,TF}, B₁ᵀ₂cols::Vector{TF}, B₂₂rows::Vector{TF}) where {TU,TF}
        _f_buf = TF()
        _f_buf2 = TF()
        Nb = length(B₁ᵀ₂cols) # number of bodies
        ibp = IBPoisson(L,R,E)
        S₀ = Matrix{Float64}(undef, Nb, Nb)
        _S⁻¹B₁ᵀ₂ = zeros(length(_f_buf), Nb)
        _S₀fact = _computeS₀fact(S₀, ibp._Sfact, B₁ᵀ₂cols, B₂₂rows, _S⁻¹B₁ᵀ₂)
        _ψ₀_buf = zeros(Nb)
        _sol_buf = IBPoissonSolution{TU,TF}(TU(),TF())
        new{Nb,TU,TF}(ibp, B₁ᵀ₂cols, B₂₂rows, _f_buf, _f_buf2, _S₀fact, _ψ₀_buf, _sol_buf)
    end
end

struct SteadyRegularizedIBPoisson{Nb,Ne,TU,TF} <: AbstractPotentialFlowSystem{TU}
    cibp::ConstrainedIBPoisson{Nb,TU,TF}
    one_vec::Vector{TF}
    e_vec::Vector{TF}
    f₀::TF
    f₀_vec::Vector{TF}

    function SteadyRegularizedIBPoisson(L::CartesianGrids.Laplacian, R::RegularizationMatrix{TU,TF}, E::InterpolationMatrix{TU,TF}, one_vec::Vector{TF}, e_vec::Vector{TF}) where {TU,TF}

        Nb = length(one_vec) # number of bodies
        Ne = length(e_vec) # number of edges

        _TU_buf = TU()

        f₀ = TF()
        _TF_buf = TF()
        _TF_buf .= 1.0
        ibp = IBPoisson(L, R, E)
        ldiv!(IBPoissonSolution(_TU_buf,f₀),ibp,IBPoissonRHS(_TU_buf,_TF_buf), onlyf=true, zerow=true)
        f₀_vec = [TF() for i=1:Nb]
        _computef₀_vec!(f₀_vec, ibp, one_vec, _TU_buf)

        R̃ = deepcopy(R)
        R̃.M .= R̃.M*Diagonal(f₀)
        cibp = ConstrainedIBPoisson(L, R̃, E, one_vec, e_vec)
        new{Nb,Ne,TU,TF}(cibp, one_vec, e_vec, f₀, f₀_vec)
    end
end

struct UnsteadyRegularizedIBPoisson{Nb,Ne,TU,TF} <: AbstractPotentialFlowSystem{TU}
    ibp::IBPoisson{TU,TF}
    one_vec::Vector{TF}
    e_vec::Vector{TF}
    vidx_vec::Vector{Vector{Int}}
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

    function UnsteadyRegularizedIBPoisson(L::CartesianGrids.Laplacian, R::RegularizationMatrix{TU,TF}, E::InterpolationMatrix{TU,TF}, one_vec::Vector{TF}, e_vec::Vector{TF}, vidx_vec::Vector{Vector{Int}}) where {TU,TF}
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
        _computef₀_vec!(f₀_vec, ibp, one_vec, _TU_buf)

        R̃ = deepcopy(R)
        R̃.M .= R̃.M*Diagonal(f₀)
        ibp = IBPoisson(L, R̃, E)

        f̃₀_vec = [TF() for i=1:Nb]
        _computef₀_vec!(f̃₀_vec, ibp, one_vec, _TU_buf)


        d_vec = [TU() for i=1:Ne]
        f̃_vec = [TF() for i=1:Ne]
        Γ₀ = sum(f₀)
        _activef̃lim_vec = zeros(Ne)
        _Souter = zeros(Nb+Ne,Nb+Ne)
        _r₂_buf = zeros(Nb+Ne)
        _y_buf = zeros(Nb+Ne)

        new{Nb,Ne,TU,TF}(ibp, one_vec, e_vec, vidx_vec, f₀, f₀_vec, f̃₀_vec, f̃_vec, Γ₀, d_vec, _activef̃lim_vec, _TU_buf, _f_buf, _Souter, _r₂_buf, _y_buf)
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

function ldiv!(sol::TS, sys::ConstrainedIBPoisson{Nb,TU,TF}, rhs::ConstrainedIBPoissonRHS) where {Nb,TU,TF,TS}

    # fstar = S⁻¹(ψb+EL⁻¹)
    ldiv!(sol, sys.ibp, rhs, onlyf=true)

    # ψ₀ = S₀\(Γb-B₂₂*fstar)
    for b in 1:Nb
        sys._ψ₀_buf[b] = dot(sys.B₂₂rows[b]',sol.f)
    end
    sys._ψ₀_buf .= rhs.fconstraintRHS .- sys._ψ₀_buf
    ldiv!(sol.ψ₀, sys._S₀fact, sys._ψ₀_buf)

    # f = fstar - S⁻¹(B₁ᵀ₂*ψ₀)
    sys._f_buf .= 0.0
    _addlincomboto!(sys._f_buf, sol.ψ₀, sys.B₁ᵀ₂cols)
    ldiv!(sys._f_buf2, sys.ibp._Sfact, sys._f_buf)
    sol.f .-= sys._f_buf2

    # ψ = L⁻¹(-w - Rf)
    ldiv!(sol, sys.ibp, rhs, useprovidedf=true)

end

function ldiv!(sol::TS, sys::SteadyRegularizedIBPoisson{Nb,Ne,TU,TF}, rhs::TR) where {Nb,Ne,TU,TF,TS,TR}

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

function _computeS₀fact(S₀,Sfact,B₁ᵀ₂cols,B₂₁rows,_S⁻¹B₁ᵀ₂)
    for c in 1:length(B₁ᵀ₂cols)
        ldiv!(view(_S⁻¹B₁ᵀ₂,:,c),Sfact,B₁ᵀ₂cols[c])
    end

    for r in 1:length(B₂₁rows)
        mul!(view(S₀,r:r,:),B₂₁rows[r]',_S⁻¹B₁ᵀ₂) # need to find better solution for this indexing
    end

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
        for vid in sys.vidx_vec[i]
            sys._Souter[end-Nb+i,end-Ne+vid] += 1.0
        end
    end

    ldiv!(sys._y_buf, factorize(sys._Souter), sys._r₂_buf)
    sol.ψ₀ .= sys._y_buf[1:Nb]
    sol.δΓ_vec .= sys._y_buf[Nb+1:end]
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

function _computef₀_vec!(f₀_vec::Vector{TF}, ibp::IBPoisson{TU,TF}, one_vec::Vector{TF}, TU_buf::TU) where {TU,TF}
    for b in 1:length(f₀_vec)
        ldiv!(IBPoissonSolution(TU_buf,f₀_vec[b]),ibp,IBPoissonRHS(TU_buf,one_vec[b]), onlyf=true, zerow=true)
    end
end
