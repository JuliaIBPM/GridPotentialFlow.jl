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

include("edges.jl")
include("bodies.jl")

# include("solver/problem.jl")
# include("solver/righthandside.jl")
# include("solver/solution.jl")
# include("solver/systems.jl")

include("vortexelements.jl")
# include("vortexmodel.jl")
# include("pressure.jl")

@ilmproblem PotentialFlow vector

struct PotentialFlowCache{HLMT,VORT,SNKT,EDGT,ST,FT,GT,VT} <: ImmersedLayers.AbstractExtraILMCache
    # Cache for Helmholtz decomposition
    helmcache :: HLMT
    # Cache for vortex elements
    vortcache :: VORT
    # Cache for sink/source elements
    sinkcache :: SNKT
    # Cache for regularizing sharp edges
    edgecache :: EDGT
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
    edgecache = nothing
    RTLinvR = create_RTLinvR_scalar(base_cache)
    γtemp = zeros_surfacescalar(base_cache)
    σtemp = zeros_surfacescalar(base_cache)
    vtemp = zeros_grid(base_cache)

    PotentialFlowCache(helmcache,vortcache,sinkcache,edgecache,RTLinvR,γtemp,σtemp,vtemp)
end

function solve!(ϕ::Nodes{Primal},ϕcache::ScalarPotentialCache,sys::ILMSystem) where {SNKT}
    println("solve ϕ system")
end

function solve!(ψ::Nodes{Dual},ψcache::VectorPotentialCache,edgecache::EDGT,sys::ILMSystem) where {EDGT<:EdgeCache}
    println("solve ψ system with sharp edges")
end

function solve!(ψ::Nodes{Dual},ψcache::VectorPotentialCache,edgecache::EDGT,sys::ILMSystem) where {EDGT<:Nothing}
    println("solve ψ system without sharp edges")

    @unpack extra_cache, base_cache, bc = sys
    @unpack RTLinvR, γtemp, vtemp = extra_cache
    @unpack dvn, vn, stemp, curlv_temp = ψcache

    prescribed_surface_jump!(dvn,sys)
    prescribed_surface_average!(vn,sys)

    regularize_normal_cross!(vtemp,dvn,base_cache)
    curl!(curlv_temp,vtemp,base_cache)
    # apply_forcing!(curlv_temp,ψ,t,fcache,phys_params)
    vectorpotential_from_curlv!(stemp,curlv_temp,base_cache)

    interpolate!(γtemp,stemp,base_cache)
    γtemp .= vn .- γtemp
    γtemp .= RTLinvR\γtemp # Should use ldiv!(γtemp,RTLinvR) for no allocation or even better implement CG so we don't have to construct RTLinvR.

    regularize!(curlv_temp,γtemp,base_cache)
    inverse_laplacian!(ψ,curlv_temp,base_cache)
    ψ .= stemp .- ψ

    return ψ
end

end
