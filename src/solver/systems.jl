import LinearAlgebra: \, ldiv!, mul!, norm, LU, factorize, Factorization, dot

using LinearMaps

abstract type AbstractPotentialFlowSystem{TU} end

"""
$(TYPEDEF)

Defines an immersed boundary Poisson system that solves for the streamfunction and bound vortex sheet strength, given a streamfunction boundary condition on immersed boundary points.

# Fields

$(TYPEDFIELDS)
"""
struct IBPoisson{TU,TF} <: AbstractPotentialFlowSystem{TU}
    """L: Discrete Laplacian."""
    L::CartesianGrids.Laplacian
    """R: Regularization operator."""
    R::RegularizationMatrix{TU,TF}
    """E: Interpolation operator."""
    E::InterpolationMatrix{TU,TF}
    """Sfact: factorized version of the Schur complement."""
    Sfact::LU # Factorization

    """Buffers"""
    _A⁻¹r₁::TU
    _B₁ᵀf::TU
    _f_buf::TF
end

function IBPoisson(L::CartesianGrids.Laplacian, R::RegularizationMatrix{TU,TF}, E::InterpolationMatrix{TU,TF}) where {TU,TF}
    _A⁻¹r₁ = TU()
    _B₁ᵀf = TU()
    _f_buf = TF()

    Rmap = LinearMap(x->vec(R*TF(x)),length(_A⁻¹r₁),length(_f_buf));
    Emap = LinearMap(x->vec(E*TU(x)),length(_f_buf),length(_A⁻¹r₁))
    L⁻¹map = LinearMap(x->vec(L\TU(x)),length(_A⁻¹r₁));
    Smap = -Emap*L⁻¹map*Rmap
    Sfact = factorize(Matrix(Smap))

    IBPoisson{TU,TF}(L, R, E, Sfact, _A⁻¹r₁, _B₁ᵀf, _f_buf)
end

"""
$(TYPEDEF)

Defines an immersed boundary Poisson system with a single type of constraints on the bound vortex sheet strength.

# Fields

$(TYPEDFIELDS)
"""
struct ConstrainedIBPoisson{Nb,TU,TF} <: AbstractPotentialFlowSystem{TU}
    """ibp: Immersed boundary poisson inner system that solves for the streamfunction and bound vortex sheet strength, given a streamfunction boundary condition on immersed boundary points."""
    ibp::IBPoisson{TU,TF}
    """B₁ᵀ₂cols: Columns of the second partition of the B₁ᵀ matrix for the outer system."""
    B₁ᵀ₂cols::Vector{TF}
    """B₂₂rows: Rows of the second partition of the B₂ matrix for the outer system."""
    B₂₂rows::Vector{TF}
    """Soutfact: Factorized version of the Schur complement for the outer system."""
    Soutfact::Union{Factorization,Float64}

    """Buffers"""
    _f_buf::TF
    _ψ₀_buf::Vector{Float64}
    _sol_buf::IBPoissonSolution{TU,TF}
end

function ConstrainedIBPoisson(L::CartesianGrids.Laplacian, R::RegularizationMatrix{TU,TF}, E::InterpolationMatrix{TU,TF}, B₁ᵀ₂cols::Vector{TF}, B₂₂rows::Vector{TF}) where {TU,TF}
    _f_buf = TF()
    Nb = length(B₁ᵀ₂cols) # number of bodies
    ibp = IBPoisson(L,R,E)
    S₀ = Matrix{Float64}(undef, Nb, Nb)
    _S⁻¹B₁ᵀ₂ = zeros(length(_f_buf), Nb)
    Soutfact = _computeSoutfact(S₀, ibp.Sfact, B₁ᵀ₂cols, B₂₂rows, _S⁻¹B₁ᵀ₂)
    _ψ₀_buf = zeros(Nb)
    _sol_buf = IBPoissonSolution{TU,TF}(TU(),TF())
    ConstrainedIBPoisson{Nb,TU,TF}(ibp, B₁ᵀ₂cols, B₂₂rows, Soutfact, _f_buf, _ψ₀_buf, _sol_buf)
end

"""
$(TYPEDEF)

Defines an immersed boundary Poisson system with constraints on the bound vortex sheet strength that regularize up to one edge per body.

# Fields

$(TYPEDFIELDS)
"""
struct SteadyRegularizedIBPoisson{Nb,Ne,TU,TF} <: AbstractPotentialFlowSystem{TU}
    """cibp: Constrained immersed boundary Poisson system."""
    cibp::ConstrainedIBPoisson{Nb,TU,TF}
    """one_vec: Array of vectors with the i-th entry of the j-th vector equal to one if i is an index of the points that belong to the j-th body and zero otherwise."""
    one_vec::Vector{TF}
    """e_vec: Array of unit vectors with the inner product of the k-th unit vector and vortex sheet strength equal to the value of the vortex sheet strength at the k-th edge."""
    e_vec::Vector{TF}
    """f₀: Unregularized bound vortex sheet strength due to a uniform unit boundary condition for the streamfunction on all immersed bodies."""
    f₀::TF
    """f₀_vec: Array of vectors with the i-th entry of the j-th vector equal to f₀[i] if i is an index of the points that belong to the j-th body and zero otherwise."""
    f₀_vec::Vector{TF}
end

function SteadyRegularizedIBPoisson(L::CartesianGrids.Laplacian, R::RegularizationMatrix{TU,TF}, E::InterpolationMatrix{TU,TF}, one_vec::Vector{TF}, e_vec::Vector{TF}) where {TU,TF}

    Nb = length(one_vec) # number of bodies
    Ne = length(e_vec) # number of edges

    _TU_buf = TU()

    f₀ = TF()
    _TF_buf = TF()
    _TF_buf .= 1.0
    ibp = IBPoisson(L, R, E)
    ldiv!(IBPoissonSolution(_TU_buf,f₀),ibp,IBPoissonRHS(_TU_buf,_TF_buf), onlyf=true, zerow=true)
    f₀_vec = [f₀.*one_vec[i] for i=1:Nb]

    R̃ = deepcopy(R)
    R̃.M .= R̃.M*Diagonal(f₀)
    cibp = ConstrainedIBPoisson(L, R̃, E, one_vec, e_vec)
    SteadyRegularizedIBPoisson{Nb,Ne,TU,TF}(cibp, one_vec, e_vec, f₀, f₀_vec)
end

"""
$(TYPEDEF)

Defines an immersed boundary Poisson system that, in addition to solving for the streamfunction and bound vortex sheet strength, solves for N point vortex strengths to satisfy N constraints on the bound vortex sheet strength that regularize specified edges.

# Fields

$(TYPEDFIELDS)
"""
struct UnsteadyRegularizedIBPoisson{Nb,Ne,TU,TF} <: AbstractPotentialFlowSystem{TU}
    """ibp: Immersed boundary poisson inner system with re-scaled regularization operator that solves for the streamfunction and bound vortex sheet strength, given a streamfunction boundary condition on immersed boundary points."""
    ibp::IBPoisson{TU,TF}
    """one_vec: Array of vectors with the i-th entry of the j-th vector equal to one if i is an index of the points that belong to the j-th body and zero otherwise."""
    one_vec::Vector{TF}
    """e_vec: Array of unit vectors with the inner product of the i-th unit vector and vortex sheet strength equal to the value of the vortex sheet strength at the i-th edge."""
    e_vec::Vector{TF}
    """vidx_vec: Array of vectors with the i-th vector containing the indices of the vortices that are shedded from the i-th body."""
    vidx_vec::Vector{Vector{Int}}
    """f₀: Unregularized vortex sheet strength due to a uniform unit boundary condition for the streamfunction on all immersed bodies."""
    f₀::TF
    """f₀_vec: Array of vectors with the i-th entry of the j-th vector equal to f₀[i] if i is an index of the points that belong to the j-th body and zero otherwise."""
    f₀_vec::Vector{TF}
    """f̃₀_vec: Array of unregularized bound vortex sheet strengths, obtained using the re-scaled regularization operator, with the i-th entry the bound vortex sheet strength due to a uniform unit boundary condition for the streamfunction on the i-th immersed body and zero everywhere else."""
    f̃₀_vec::Vector{TF}
    """f̃_vec: Array of unregularized bound vortex sheet strengths, obtained using the re-scaled regularization operator, with the i-th entry the bound vortex sheet strength due to the i-th vorticity field in `d_vec` and a zero boundary condition for the streamfunction on all bodies."""
    f̃_vec::Vector{TF}
    """Γ₀: The circulation of the of the bound vortex sheet strength f₀."""
    Γ₀::Float64
    """d_vec: Array of vorticity fields with the i-th entry the vorticity field due to the i-th point vortex that is used for regularizing the edges with its strength set to one and the other point vortex strengths set to zero."""
    d_vec::Vector{TU}
    """activef̃lim_vec: Array of f̃ values that are used as the right-hand side for the constraints on f̃ when solving the system. These values are chosen by comparing the unconstrained f̃ from the inner system to the limits provided in `rhs` to `ldiv!` and selecting the limit that is exceeded or setting it to `Inf` if no limit is exceeded."""
    activef̃lim_vec::Vector{f̃Limit}
    """Sout: Schur complement for the outer system."""
    Sout::Matrix{Float64}

    """Buffers"""
    _TU_buf::TU
    _f_buf::TF
    _r₂_buf::Vector{Float64}
    _y_buf::Vector{Float64}
end

function UnsteadyRegularizedIBPoisson(L::CartesianGrids.Laplacian, R::RegularizationMatrix{TU,TF}, E::InterpolationMatrix{TU,TF}, one_vec::Vector{TF}, e_vec::Vector{TF}, vidx_vec::Vector{Vector{Int}}) where {TU,TF}
    Nb = length(one_vec) # number of bodies
    Ne = length(e_vec) # number of edges

    ibp = IBPoisson(L, R, E)

    _TU_buf = TU()
    _f_buf = TF()

    f₀ = TF()

    _f_buf .= 1.0
    ldiv!(IBPoissonSolution(_TU_buf,f₀),ibp,IBPoissonRHS(_TU_buf,_f_buf), onlyf=true, zerow=true)

    f₀_vec = [f₀.*one_vec[i] for i=1:Nb]

    R̃ = deepcopy(R)
    R̃.M .= R̃.M*Diagonal(f₀)
    ibp = IBPoisson(L, R̃, E)

    f̃₀_vec = [TF() for i=1:Nb]
    _computef₀_vec!(f̃₀_vec, ibp, one_vec, _TU_buf)


    d_vec = [TU() for i=1:Ne]
    f̃_vec = [TF() for i=1:Ne]
    Γ₀ = sum(f₀)
    activef̃lim_vec = zeros(Ne)
    Sout = zeros(Nb+Ne,Nb+Ne)
    _r₂_buf = zeros(Nb+Ne)
    _y_buf = zeros(Nb+Ne)

    UnsteadyRegularizedIBPoisson{Nb,Ne,TU,TF}(ibp, one_vec, e_vec, vidx_vec, f₀, f₀_vec, f̃₀_vec, f̃_vec, Γ₀, d_vec, activef̃lim_vec, Sout, _TU_buf, _f_buf, _r₂_buf, _y_buf)
end

function ldiv!(sol::TS, sys::IBPoisson{TU,TF}, rhs::TR; zerow=false, zeroψb=false, useprovidedf=false, onlyf=false) where {TU,TF,TS<:AbstractIBPoissonSolution,TR<:AbstractIBPoissonRHS}

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

        ldiv!(sol.f, sys.Sfact, sys._f_buf) # S⁻¹(r₂ - B₂A⁻¹r₁)

        if !onlyf # ψ = L⁻¹(-w-Rf)
            mul!(sys._B₁ᵀf, sys.R, sol.f)
            ldiv!(sol.ψ, sys.L, sys._B₁ᵀf)
            sol.ψ .= sys._A⁻¹r₁ .- sol.ψ # L⁻¹(-w) was already calculated earlier and stored in sys._A⁻¹r₁
        end
    end

    return sol
end

function ldiv!(sol::ConstrainedIBPoissonSolution{TU,TF}, sys::ConstrainedIBPoisson{Nb,TU,TF}, rhs::ConstrainedIBPoissonRHS) where {Nb,TU,TF}

    # fstar = S⁻¹(ψb+EL⁻¹w)
    ldiv!(sol, sys.ibp, rhs, onlyf=true)

    # ψ₀ = S₀\(Γb-B₂₂*fstar)
    for b in 1:Nb
        sys._ψ₀_buf[b] = dot(sys.B₂₂rows[b]',sol.f)
    end
    sys._ψ₀_buf .= rhs.fconstraintRHS .- sys._ψ₀_buf
    ldiv!(sol.ψ₀, sys.Soutfact, sys._ψ₀_buf)

    # f = fstar - S⁻¹(B₁ᵀ₂*ψ₀)
    sys._f_buf .= 0.0
    _addlincomboto!(sys._f_buf, sol.ψ₀, sys.B₁ᵀ₂cols)
    ldiv!(sys.ibp.Sfact, sys._f_buf)
    sol.f .-= sys._f_buf

    # ψ = L⁻¹(-w - Rf)
    ldiv!(sol, sys.ibp, rhs, useprovidedf=true)

end

function ldiv!(sol::ConstrainedIBPoissonSolution{TU,TF}, sys::SteadyRegularizedIBPoisson{Nb,Ne,TU,TF}, rhs::ConstrainedIBPoissonRHS) where {Nb,Ne,TU,TF}

    ldiv!(sol, sys.cibp, rhs)
    sol.f .*= sys.f₀

end

function ldiv!(sol::ConstrainedIBPoissonSolution{TU,TF}, sys::UnsteadyRegularizedIBPoisson{Nb,Ne,TU,TF}, rhs::UnsteadyRegularizedIBPoissonRHS{TU,TF}) where {Nb,Ne,TU,TF}

    # f̃star = S̃⁻¹(ψb+EL⁻¹w)
    ldiv!(sol, sys.ibp, rhs, onlyf=true)

    for k in 1:Ne
        # eq 2.47
        ldiv!(IBPoissonSolution(sys._TU_buf, sys.f̃_vec[k]), sys.ibp, IBPoissonRHS(sys.d_vec[k],sys._f_buf), zeroψb=true, onlyf=true)
        # eq 2.60
        sys.activef̃lim_vec[k] = _findactivef̃limit(sys.e_vec[k], sol.f, rhs.f̃lim_vec[k])
    end

    # The shedding edges are the ones for which the active f̃ limit is not Inf
    sheddingedges = [e for e in 1:length(sys.activef̃lim_vec) if sys.activef̃lim_vec[e] != Inf]

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

function _computeSoutfact(S₀,Sfact,B₁ᵀ₂cols,B₂₁rows,_S⁻¹B₁ᵀ₂)
    for c in 1:length(B₁ᵀ₂cols)
        ldiv!(view(_S⁻¹B₁ᵀ₂,:,c),Sfact,B₁ᵀ₂cols[c])
    end

    for r in 1:length(B₂₁rows)
        mul!(view(S₀,r:r,:),B₂₁rows[r]',_S⁻¹B₁ᵀ₂) # need to find better solution for this indexing
    end

    S₀ .= .-S₀

    Soutfact = factorize(S₀)
end

function ldiv!(x::ScalarData ,A::LU, f::ScalarData)
    ldiv!(x.data,A,f.data)
end

function ldiv!(A::LU, f::ScalarData)
    ldiv!(A,f.data)
end

"""
$(SIGNATURES)

Computes the vortex strengths for the edges `sheddingedges` and ψ₀ values and stores them in the `δΓ_vec` and `ψ₀` fields of `sol`. If there are edges that are not in `sheddingedges`, their corresponding vortex strengths are set to zero.
"""
function _computeδΓandψ₀!(sol::ConstrainedIBPoissonSolution{TU,TF}, sys::UnsteadyRegularizedIBPoisson{Nb,Ne,TU,TF}, rhs::UnsteadyRegularizedIBPoissonRHS{TU,TF}, sheddingedges::Vector{Int}) where {Nb,Ne,TU,TF}
    sys.Sout .= 0.0
    for i in 1:Ne
        if i in sheddingedges
            for b in 1:Nb
                sys.Sout[i,b] = dot(sys.e_vec[i], sys.f̃₀_vec[b])
            end
            for e in 1:Ne
                sys.Sout[i,end-Ne+e] = -dot(sys.e_vec[i], sys.f̃_vec[e])
            end
            sys._r₂_buf[i] = sys.activef̃lim_vec[i] - dot(sys.e_vec[i], sol.f)
        else
            sys.Sout[i,end-Ne+i] = 1.0
            sys._r₂_buf[i] = 0.0
        end
    end
    for i in 1:Nb
        for b in 1:Nb
            sys.Sout[end-Nb+i,b] = dot(sys.f₀_vec[i], sys.f̃₀_vec[b])
        end
        for e in 1:Ne
            sys.Sout[end-Nb+i,end-Ne+e] = -dot(sys.f₀_vec[i], sys.f̃_vec[e])
        end
        sys._r₂_buf[end-Nb+i] = -rhs.Γw_vec[i] - dot(sys.f₀_vec[i], sol.f)
    end
    sys.Sout .= .-sys.Sout
    for i in 1:Nb
        for vid in sys.vidx_vec[i]
            sys.Sout[end-Nb+i,end-Ne+vid] += 1.0
        end
    end

    ldiv!(sys._y_buf, factorize(sys.Sout), sys._r₂_buf)
    sol.ψ₀ .= sys._y_buf[1:Nb]
    sol.δΓ_vec .= sys._y_buf[Nb+1:end]
end

function _findactivef̃limit(e::AbstractVector, f̃::TF, f̃lim::f̃Limits) where {TF}

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
