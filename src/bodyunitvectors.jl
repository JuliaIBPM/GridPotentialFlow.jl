struct BodyUnitVector{TF,TK}
    data::TF
    k::TK
    closuretype::Type{<:RigidBodyTools.BodyClosureType}
end

function BodyUnitVector(Nk::Integer,k::Integer,closuretype::Type{<:RigidBodyTools.BodyClosureType})

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

    return BodyUnitVector{typeof(data),typeof(k)}(data,k,closuretype)
end
