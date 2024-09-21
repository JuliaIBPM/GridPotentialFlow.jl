struct GridPotentialILMProblem{DT,ST,DTP,PHT,BCF,FF,DTF,MTF} <: AbstractScalarILMProblem{DT,ST,DTP}
   g :: PhysicalGrid
   bodies :: BodyList
   phys_params :: PHT
   bc :: BCF
   forcing :: FF
   timestep_func :: DTF
   motions :: MTF
   GridPotentialILMProblem(g::PT,bodies::BodyList;ddftype=CartesianGrids.Yang3,scaling=IndexScaling,reference_body=0,phys_params=nothing,bc=nothing,forcing=nothing,timestep_func=nothing,motions=nothing) where {PT} = (motion_data = ImmersedLayers._construct_motion(motions,reference_body,bodies); new{ddftype,scaling,ImmersedLayers.DEFAULT_DATA_TYPE,typeof(phys_params),typeof(bc),typeof(forcing),typeof(timestep_func),typeof(motion_data)}(g,bodies,phys_params,bc,forcing,timestep_func,motion_data))
   GridPotentialILMProblem(g::PT,body::Body;ddftype=CartesianGrids.Yang3,scaling=IndexScaling,reference_body=0,phys_params=nothing,bc=nothing,forcing=nothing,timestep_func=nothing,motions=nothing) where {PT} = (motion_data = ImmersedLayers._construct_motion(motions,reference_body,BodyList([body])); new{ddftype,scaling,ImmersedLayers.DEFAULT_DATA_TYPE,typeof(phys_params),typeof(bc),typeof(forcing),typeof(timestep_func),typeof(motion_data)}(g,BodyList([body]),phys_params,bc,forcing,timestep_func,motion_data))
end

GridPotentialILMProblem(g,body::PotentialFlowBody;kwargs...) = GridPotentialILMProblem(g,body.points;kwargs...)
function GridPotentialILMProblem(g,bodies::Vector{<:PotentialFlowBody};kwargs...)
    bl = BodyList()
    for b in bodies
        push!(bl,b.points)
    end
    GridPotentialILMProblem(g,bl;kwargs...)
end


struct GridPotentialExtraCache{CT,RNT,ENT,RCT,ECT,VT,CVT} <: AbstractExtraILMCache
   CLinvCT :: CT
   Rn :: RNT
   En :: ENT
   Rc :: RCT
   Ec :: ECT
   v_cache :: VT
   cv_cache :: CVT
end

function ImmersedLayers.prob_cache(prob::GridPotentialILMProblem,base_cache::BasicILMCache)
    @unpack g, regop,gcurl_cache,gdata_cache,sdata_cache  = base_cache
    A = ImmersedLayers.create_CLinvCT(base_cache,scale=cellsize(g))
    Rn = RegularizationMatrix(regop,sdata_cache,gcurl_cache)
    En = InterpolationMatrix(regop,gcurl_cache,sdata_cache)
    Rc = RegularizationMatrix(regop,sdata_cache,gdata_cache)
    Ec = InterpolationMatrix(regop,gdata_cache,sdata_cache)
    v_cache = Edges(Primal,size(g))
    cv_cache = ConvectiveDerivativeCache(EdgeGradient(Primal,size(g)))

    GridPotentialExtraCache(A,Rn,En,Rc,Ec,v_cache,cv_cache)
end
