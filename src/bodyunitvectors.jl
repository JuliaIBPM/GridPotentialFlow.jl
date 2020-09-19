import SparseArrays: spzeros

export BodyUnitVector

const BodyUnitVector = ScalarData

function BodyUnitVector(Nk::Integer,k::Integer,closuretype::Type{<:RigidBodyTools.BodyClosureType};midpoints::Bool=false)

    if midpoints
        @assert 1 <= k <= Nk+1 "k has to be in the range 1:Nk+1"
    else
        @assert 1 <= k <= Nk "k has to be in the range 1:Nk"
    end

    data = spzeros(Nk)

    # Change to something based on a coordinate?
    if midpoints
        if closuretype == RigidBodyTools.OpenBody
            if k == 1
                data[1] = 1.5
                data[2] = -0.5
            elseif k == Nk+1
                data[Nk-1] = -0.5
                data[Nk] = 1.5
            else
                data[k-1] = 0.5 # I'm not sure if these values are correct
                data[k] = 0.5
            end
        else
            error
        end
    else
        data[k] = 1.0;
    end

    return BodyUnitVector(data)
end

function BodyUnitVector(body::Body{N,C},k::Integer;midpoints::Bool=false) where {N,C}
    return BodyUnitVector(length(body),k,C,midpoints=midpoints)
end
