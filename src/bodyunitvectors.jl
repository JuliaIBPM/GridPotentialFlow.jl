struct BodyUnitVector{Nk,T,Tk}
    data::ScalarData{Nk,T}
    k::Tk
    closuretype::Type{<:RigidBodyTools.BodyClosureType}
end

function BodyUnitVector(Nk::Int,k::Int,closuretype::Type{<:RigidBodyTools.BodyClosureType})

    @assert 1 <= k <= Nk "k has to be in the range 1:Nk"

    data = ScalarData(Nk)

    # ONLY VALID IF POINTS ARE MIDPOINTS
    if closuretype == RigidBodyTools.OpenBody
        if k == 1
            data[1] = 1.5
            data[2] = -0.5
        elseif k == Nk
            data[Nk-1] = -0.5
            data[Nk] = 1.5
        else
            data[k-1] = 0.5
            data[k] = 0.5
        end
    else
        error
    end

    return BodyUnitVector{Nk,Float64,Int64}(data,k,closuretype)
end
