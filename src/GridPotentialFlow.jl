module GridPotentialFlow

using Reexport
using StructArrays
using DocStringExtensions
using RecipesBase
using UnPack
using LinearAlgebra
import MacroTools.@forward

@reexport using ImmersedLayers

export PotentialFlowProblem, setup_problem, streamfunction, scalarpotential, impulse

include("bodies.jl")
include("constraints.jl")
include("vortexelements.jl")
# include("vortexmodel.jl")
# include("pressure.jl")

#= Supporting functions =#

setup_problem(g;kwargs...) =
    PotentialFlowProblem(g;kwargs...)

setup_problem(g,bl;bc=nothing,kwargs...) =
    PotentialFlowProblem(g,bl;bc=get_bc_func(bc),kwargs...)

function default_freestream(t,phys_params)
    Vinfmag = get(phys_params,"freestream speed",0.0)
    Vinf_angle = get(phys_params,"freestream angle",0.0)
    Uinf = Vinfmag*cos(Vinf_angle)
    Vinf = Vinfmag*sin(Vinf_angle)
    return Uinf, Vinf
end

function default_vbplus(t,base_cache,phys_params,motions)
  vbplus = zeros_surface(base_cache)
  return vbplus
end

function default_vbminus(t,base_cache,phys_params,motions)
    vbminus = zeros_surface(base_cache)
    return vbminus
end

const DEFAULT_DS_TO_DX_RATIO = 1.4
const DEFAULT_FREESTREAM_FUNC = default_freestream
const DEFAULT_VBPLUS_FUNC = default_vbplus
const DEFAULT_VBMINUS_FUNC = default_vbminus

#=
Process keywords
=#

function get_freestream_func(forcing::Dict)
    return get(forcing,"freestream",DEFAULT_FREESTREAM_FUNC)
end

get_freestream_func(::Nothing) = get_freestream_func(Dict())

function get_forcing_models(forcing::Dict)
    return get(forcing,"forcing models",nothing)
end

get_forcing_models(::Nothing) = get_forcing_models(Dict())

function get_bc_func(bc_in::Dict)
    bc = Dict()
    bc["exterior"] = haskey(bc_in,"exterior") ? bc_in["exterior"] : DEFAULT_VBPLUS_FUNC
    bc["interior"] = haskey(bc_in,"interior") ? bc_in["interior"] : DEFAULT_VBMINUS_FUNC
    return bc
end

get_bc_func(::Nothing) = get_bc_func(Dict())

#=
Defining the extra cache and extending prob_cache
=#

@ilmproblem PotentialFlow vector

struct PotentialFlowcache{HLMT,VORT,SNKT,CONT,PHIST,PSIST,FT,GT,VT,VBT} <: ImmersedLayers.AbstractExtraILMCache
    # Cache for Helmholtz decomposition
    helmcache :: HLMT
    # Cache for vortex elements
    vortcache :: VORT
    # Cache for sink/source elements
    sinkcache :: SNKT
    # Cache for regularizing sharp edges
    constraintscache :: CONT
    # CLinvCT
    CLinvCT :: PHIST
    # RTLinvR
    RTLinvR :: PSIST
    # Vortex sheet strength
    γtemp :: FT
    # Jump in scalar potential on the body
    dϕtemp :: GT
    # Velocity
    vtemp :: VT
    # Body velocity
    vbtemp :: VBT
end

function ImmersedLayers.prob_cache(prob::PotentialFlowProblem,base_cache::BasicILMCache{N,scaling}) where {N,scaling}
    @unpack bl = base_cache

    helmcache = VectorFieldCache(base_cache)
    vortcache = nothing
    sinkcache = nothing
    CLinvCT = create_CLinvCT_scalar(base_cache)
    RTLinvR = create_RTLinvR_scalar(base_cache)
    γtemp = zeros_surfacescalar(base_cache)
    dϕtemp = zeros_surfacescalar(base_cache)
    vtemp = zeros_grid(base_cache)
    vbtemp = zeros_surface(base_cache)

    constraintscache = ConstraintsCache(RTLinvR,base_cache)

    PotentialFlowcache(helmcache,vortcache,sinkcache,constraintscache,CLinvCT,RTLinvR,γtemp,dϕtemp,vtemp,vbtemp)
end

function solve!(ψ::Nodes{Dual},γ,ψb,dψb,constraintscache::CONT,sys::ILMSystem,t) where {CONT<:ConstraintsCache}
    println("solve non-shedding ψ system with extra constraints")

    @unpack extra_cache, base_cache, bc = sys
    @unpack helmcache, constraintscache, RTLinvR, vtemp = extra_cache
    @unpack stemp, curlv_temp = helmcache.wcache
    @unpack ψ₀_vec, f₀_vec, B₂₂_vec, S = constraintscache

    Nb = length(base_cache.bl)

    # x⋆ = A⁻¹r₁
    solve!(ψ,γ,ψb,dψb,nothing,sys,t)
    # y = S⁻¹(r₂-B₂x⋆)
    rhs = zeros(Nb)
    for i in 1:Nb
        rhs[i] = constraint(base_cache.bl[i]) - B₂₂_vec[i]'*γ
    end
    ψb₀ = S\rhs
    println(ψb₀)

    # x = x⋆-A⁻¹B₁ᵀy
    _subtractlincombo!(ψ,ψb₀,ψ₀_vec)
    _subtractlincombo!(γ,ψb₀,f₀_vec)
end

function solve!(ψ::Nodes{Dual},γ,ψb,dψb,constraintscache::CONT,sys::ILMSystem,t) where {CONT<:Nothing}
    println("solve ψ system without extra constraints")

    @unpack extra_cache, forcing, base_cache, phys_params = sys
    @unpack bl = base_cache
    @unpack helmcache, RTLinvR, vtemp = extra_cache
    @unpack stemp, curlv_temp = helmcache.wcache

    regularize_normal_cross!(vtemp,dψb,base_cache)
    curl!(curlv_temp,vtemp,base_cache)
    # apply_forcing!(curlv_temp,ψ,t,fcache,phys_params)
    vectorpotential_from_curlv!(stemp,curlv_temp,base_cache)

    interpolate!(γ,stemp,base_cache)
    γ .= ψb .- γ
    γ .= RTLinvR\γ # Should use ldiv!(γtemp,RTLinvR) for no allocation or even better implement CG so we don't have to construct RTLinvR.

    regularize!(curlv_temp,γ,base_cache)
    inverse_laplacian!(ψ,curlv_temp,base_cache)
    ψ .= stemp .- ψ

    return ψ
end

function solve!(ψ::Nodes{Dual},dψ,dϕ,dvn,sys::ILMSystem,t)
    @unpack extra_cache, base_cache, bc = sys
    @unpack helmcache, constraintscache, CLinvCT, vtemp = extra_cache
    @unpack ftemp, divv_temp = helmcache.dcache
    @unpack stemp, curlv_temp = helmcache.wcache

    regularize!(ftemp,dvn,base_cache)
    inverse_laplacian!(ftemp,base_cache)

    surface_curl!(stemp,dϕ,base_cache)
    surface_grad_cross!(dψ,ftemp,base_cache)
    dψ .= CLinvCT\dψ

    surface_curl_cross!(ψ,dψ,base_cache)
    ψ .-= stemp
    ψ .*= -1.0

    inverse_laplacian!(ψ,base_cache)

    return ψ
end

function solve!(ϕ::Nodes{Primal},dϕ,vn,dvn,sys::ILMSystem,t)
    @unpack extra_cache, base_cache, bc = sys
    @unpack helmcache, constraintscache, CLinvCT, vtemp = extra_cache
    @unpack ftemp, divv_temp = helmcache.dcache

    # Find the potential
    regularize!(ftemp,dvn,base_cache)
    inverse_laplacian!(ftemp,base_cache)

    surface_grad!(dϕ,ftemp,base_cache)
    dϕ .= vn - dϕ
    dϕ .= -(CLinvCT\dϕ);

    surface_divergence!(ϕ,dϕ,base_cache)
    inverse_laplacian!(ϕ,base_cache)
    ϕ .+= ftemp

    return ϕ, dϕ
end

function streamfunction(sys::ILMSystem,t)
    @unpack base_cache, extra_cache, forcing = sys
    @unpack helmcache, constraintscache, γtemp, vbtemp = extra_cache
    @unpack dv, dvn, vn, stemp = helmcache.wcache

    ψ = zeros_gridcurl(sys)
    pts = points(sys)

    # get Dirichlet boundary conditions for ψ
    surface_velocity!(vbtemp,sys,t)
    prescribed_surface_average!(dv,t,sys)
    dv .+= vbtemp
    _ψbfromvb!(vn,dv,sys)
    prescribed_surface_jump!(dv,t,sys)
    _ψbfromvb!(dvn,dv,sys)

    # subtract influence of free stream
    freestream_func = GridPotentialFlow.get_freestream_func(forcing)
    Uinf, Vinf = freestream_func(t,sys.phys_params)
    vn .-= Uinf .* pts.v .- Vinf .* pts.u

    # compute streamfunction
    solve!(ψ,γtemp,vn,dvn,constraintscache,sys,t)

    # add free stream streamfunction field
    vectorpotential_uniformvecfield!(stemp,Uinf,Vinf,base_cache)
    ψ .+= stemp

    return ψ
end

streamfunction(sys::ILMSystem) = streamfunction(sys,0.0)

function scalarpotential(sys::ILMSystem,t)
    @unpack base_cache, extra_cache, forcing = sys
    @unpack nrm = base_cache
    @unpack helmcache, constraintscache, dϕtemp, vbtemp = extra_cache
    @unpack dv, dvn, vn, ftemp = helmcache.dcache

    ϕ = zeros_griddiv(sys)
    pts = points(sys)

    # get Neumann boundary conditions for ϕ
    surface_velocity!(vbtemp,sys,t)
    prescribed_surface_average!(dv,t,sys)
    dv .+= vbtemp
    vn .= (nrm.u .* dv.u .+ nrm.v .* dv.v)
    prescribed_surface_jump!(dv,t,sys)
    dvn .= (nrm.u .* dv.u .+ nrm.v .* dv.v)

    # subtract influence of free stream
    freestream_func = GridPotentialFlow.get_freestream_func(forcing)
    Uinf, Vinf = freestream_func(t,sys.phys_params)
    vn .-= (Uinf .* nrm.u .+ Vinf .* nrm.v)

    # compute scalar potential
    solve!(ϕ,dϕtemp,vn,dvn,sys,t)

    # add free stream scalar potential field
    scalarpotential_uniformvecfield!(ftemp,Uinf,Vinf,base_cache)
    ϕ .+= ftemp
end

scalarpotential(sys::ILMSystem) = scalarpotential(sys,0.0)

function impulse(ϕ::Nodes{Primal},dϕ::ScalarData,sys,t)
    @unpack ds, nrm = sys.base_cache

    ϕb = zeros_surfacescalar(sys)

    # compute ϕ̄b
    interpolate!(ϕb,ϕ,sys)

    # compute ϕ⁺b = ϕ̄b .+ 1/2*dϕ
    ϕb .+= 1/2*dϕ

    ϕ⁺bn = zeros_surface(sys)
    ϕ⁺bn.u .= ϕb.*nrm.u
    ϕ⁺bn.v .= ϕb.*nrm.v

    P = -integrate(ϕ⁺bn,ds)

    return P
end

"""
    impulse(γ::ScalarData,sys,t)

Computes the linear impulse of the exterior flow associated with the bound vortex sheet strength `γ`, calculated as ``P = \\int_{\\mathcal{S}_\\mathrm{b}} x \\times (n \\times v^{+}) \\mathrm{d}S = \\int_{\\mathcal{S}_\\mathrm{b}} x \\times (\\gamma + n \\times v_\\mathrm{b}) \\mathrm{d}S``, with ``\\gamma = n \\times (v^{+} - v^{-})``. Note that if ``n \\times v^{-} \\neq n \\times v_\\mathrm{b}``, this calculation will be incorrect.
"""
function impulse(γ::ScalarData,sys,t)
    @unpack bl, ds, nrm, sdata_cache, sscalar_cache = sys.base_cache
    @unpack dv = sys.extra_cache.helmcache.wcache

    surface_velocity!(dv,sys,t)

    P = [0.0,0.0]
    for i in 1:length(bl)
        r = getrange(bl,i)
        P += _impulsesurfaceintegral(bl[i], γ[r], dv.u[r], dv.v[r], ds[r], nrm.u[r], nrm.v[r], sdata_cache.u[r], sdata_cache.v[r], sscalar_cache[r])
    end
    return P
end

function impulse(sys,t)

end

function _impulsesurfaceintegral(body::Body{N,RigidBodyTools.ClosedBody}, f, u, v, ds, nx, ny, xcrossncrossv_x_cache, xcrossncrossv_y_cache, ncrossv_cache) where {N}
    rx,ry = collect(body)
    ncrossv_cache .= f .+ nx .* v .- ny .* u
    xcrossncrossv_x_cache .= ncrossv_cache .* ry
    xcrossncrossv_y_cache .= .-ncrossv_cache .* rx
    return integrate(VectorData(xcrossncrossv_x_cache,xcrossncrossv_y_cache),ScalarData(ds))
end

function _impulsesurfaceintegral(body::Body{N,RigidBodyTools.OpenBody}, f, u, v, ds, nx, ny, xcrossncrossv_x_cache, xcrossncrossv_y_cache, ncrossv_cache) where {N}
    rx,ry = collect(body)
    ncrossv_cache .= f
    xcrossncrossv_x_cache .= ncrossv_cache .* ry
    xcrossncrossv_y_cache .= .-ncrossv_cache .* rx
    return integrate(VectorData(xcrossncrossv_x_cache,xcrossncrossv_y_cache),ScalarData(ds))
end

function _addlincombo!(s,c,A)
    for i in 1:length(c)
        s .+= c[i].*A[i]
    end
end

function _subtractlincombo!(s,c,A)
    for i in 1:length(c)
        s .-= c[i].*A[i]
    end
end

function _ψbfromvb!(ψb,vb,sys)
    pts = points(sys)
    ψb .= vb.u*pts.v - vb.v*pts.u
end

end
