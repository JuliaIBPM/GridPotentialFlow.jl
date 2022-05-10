module GridPotentialFlow

using Reexport
using StructArrays
using DocStringExtensions
using RecipesBase
using UnPack
using LinearAlgebra
import MacroTools.@forward

@reexport using ImmersedLayers

export PotentialFlowProblem, solve!

include("bodies.jl")
include("constraints.jl")

# include("solver/problem.jl")
# include("solver/righthandside.jl")
# include("solver/solution.jl")
# include("solver/systems.jl")

include("vortexelements.jl")
# include("vortexmodel.jl")
# include("pressure.jl")

#= Supporting functions =#

setup_problem(g;kwargs...) =
    ViscousIncompressibleFlowProblem(g;kwargs...)

setup_problem(g,bl;bc=nothing,kwargs...) =
    ViscousIncompressibleFlowProblem(g,bl;bc=get_bc_func(bc),kwargs...)

function default_freestream(t,phys_params)
    Vinfmag = get(phys_params,"freestream speed",0.0)
    Vinf_angle = get(phys_params,"freestream angle",0.0)
    Uinf = Vinfmag*cos(Vinf_angle)
    Vinf = Vinfmag*sin(Vinf_angle)
    return Uinf, Vinf
end

function default_vnplus(t,base_cache,phys_params,motions)
  vnplus = zeros_surface(base_cache)
  return vnplus
end

function default_vnminus(t,base_cache,phys_params,motions)
    vnminus = zeros_surface(base_cache)
    return vnminus
end

const DEFAULT_DS_TO_DX_RATIO = 1.4
const DEFAULT_FREESTREAM_FUNC = default_freestream
const DEFAULT_VNPLUS_FUNC = default_vnplus
const DEFAULT_VNMINUS_FUNC = default_vnminus

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
    bc["exterior"] = haskey(bc_in,"exterior") ? bc_in["exterior"] : DEFAULT_VNPLUS_FUNC
    bc["interior"] = haskey(bc_in,"interior") ? bc_in["interior"] : DEFAULT_VNMINUS_FUNC
    return bc
end

get_bc_func(::Nothing) = get_bc_func(Dict())

#=
Defining the extra cache and extending prob_cache
=#

@ilmproblem PotentialFlow vector

struct PotentialFlowCache{HLMT,VORT,SNKT,CONT,ST,FT,GT,VT} <: ImmersedLayers.AbstractExtraILMCache
    # Cache for Helmholtz decomposition
    helmcache :: HLMT
    # Cache for vortex elements
    vortcache :: VORT
    # Cache for sink/source elements
    sinkcache :: SNKT
    # Cache for regularizing sharp edges
    constraintscache :: CONT
    # RTLinvR
    RTLinvR :: ST
    # Vortex sheet strength
    γtemp :: FT
    # Sink/source sheet strength
    σtemp :: GT
    # Velocity
    vtemp :: VT
end

function ImmersedLayers.prob_cache(prob::PotentialFlowProblem,base_cache::BasicILMCache{N,scaling}) where {N,scaling}
    @unpack bl = base_cache

    helmcache = VectorFieldCache(base_cache)
    vortcache = nothing
    sinkcache = nothing
    RTLinvR = create_RTLinvR_scalar(base_cache)
    γtemp = zeros_surfacescalar(base_cache)
    σtemp = zeros_surfacescalar(base_cache)
    vtemp = zeros_grid(base_cache)

    constraintscache = ConstraintsCache(RTLinvR,base_cache)

    PotentialFlowCache(helmcache,vortcache,sinkcache,constraintscache,RTLinvR,γtemp,σtemp,vtemp)
end

function solve!(ϕ::Nodes{Primal},σ::ScalarData,ϕcache::ScalarPotentialCache,sys::ILMSystem) where {SNKT}
    println("solve ϕ system")
end

function solve!(ψ::Nodes{Dual},γ::ScalarData,ψcache::VectorPotentialCache,constraintscache::CONT,sys::ILMSystem,t) where {CONT<:ConstraintsCache}
    println("solve non-shedding ψ system with constraints")

    @unpack extra_cache, base_cache, bc = sys
    @unpack constraintscache, RTLinvR, vtemp = extra_cache
    @unpack dvn, vn, stemp, curlv_temp = ψcache
    @unpack ψ₀_vec, f₀_vec, B₂₂_vec, S = constraintscache

    Nb = length(base_cache.bl)

    # x⋆ = A⁻¹r₁
    solve!(ψ,γ,ψcache,nothing,sys,t)
    # S = -C-B₂A⁻¹B₁ᵀ
    rhs = zeros(Nb)
    for i in 1:Nb
        rhs[i] = constraint(base_cache.bl[i]) - B₂₂_vec[i]'*γ
    end
    # y = S⁻¹(r₂-B₂x⋆)
    ψb₀ = S\rhs
    println(ψb₀)

    # x = x⋆-A⁻¹B₁ᵀy
    _subtractlincombo!(ψ,ψb₀,ψ₀_vec)
    _subtractlincombo!(γ,ψb₀,f₀_vec)
end

function solve!(ψ::Nodes{Dual},γ::ScalarData,ψcache::VectorPotentialCache,constraintscache::CONT,sys::ILMSystem,t) where {CONT<:Nothing}
    println("solve ψ system without constraints")

    @unpack extra_cache, forcing, base_cache, bc, phys_params = sys
    @unpack bl = base_cache
    @unpack RTLinvR, vtemp = extra_cache
    @unpack dv, stemp, curlv_temp = ψcache

    pts = points(sys)

    prescribed_surface_jump!(dv,sys) # should use cached variables here
    dψb = dv.u*pts.v - dv.v*pts.u
    prescribed_surface_average!(dv,sys)
    ψb = dv.u*pts.v - dv.v*pts.u

    # subtract influence of free stream
    freestream_func = get_freestream_func(forcing)
    Uinf, Vinf = freestream_func(t,phys_params)
    ψb .-= Uinf .* pts.v .- Vinf .* pts.u

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

end
