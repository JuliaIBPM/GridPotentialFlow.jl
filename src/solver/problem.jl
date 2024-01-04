struct GridPotentialILMProblem{DT,ST} <: AbstractScalarILMProblem{DT,ST}
   g :: PhysicalGrid
   bodies :: BodyList
   GridPotentialILMProblem(g::PT,bodies::BodyList;ddftype=CartesianGrids.Yang3,scaling=IndexScaling) where {PT} = new{ddftype,scaling}(g,bodies)
   GridPotentialILMProblem(g::PT,body::Body;ddftype=CartesianGrids.Yang3,scaling=IndexScaling) where {PT} = new{ddftype,scaling}(g,BodyList([body]))
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
