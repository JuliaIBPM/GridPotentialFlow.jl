const SuctionParameter = Float64

struct SuctionParameterRange
    min::SuctionParameter
    max::SuctionParameter

    function SuctionParameterRange(v1,v2)
        if v1 ≤ v2
            min = v1
            max = v2
        else
            min = v2
            max = v1
        end
        new(min,max)
    end
end
