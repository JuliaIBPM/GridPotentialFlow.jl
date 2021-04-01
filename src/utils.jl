import Base: *, /

function _computeregularizationmatrix(g::PhysicalGrid, X::VectorData{N}, f::ScalarData{N}, s::Nodes; ddftype=CartesianGrids.Yang3) where {N}

    regop = Regularize(X, cellsize(g), I0=origin(g), ddftype = ddftype, issymmetric=true)
    Rmat,_ = RegularizationMatrix(regop, f, s)
    Emat = InterpolationMatrix(regop, s, f)

    return Rmat, Emat
end

function _computecrossproductsurfaceintegral(body::Body, z)
    rx,ry = collect(body)
    @assert length(rx) == length(z)
    return [ry'*z, -rx'*z]
end

function _computebodypointsvelocity!(Ubvec, Ub, bodies)

    # Ensure that Ub is an array whose length equals to the number of bodies
    if Ub isa Tuple
        Ub = [Ub]
    end
    if isnothing(Ub)
        Ub = fill((0.0,0.0),length(bodies))
    end

    @assert length(Ub) == length(bodies)

    Ubvec.u .= 0
    Ubvec.v .= 0

    for i in 1:length(bodies)
        Ubvec.u[getrange(bodies,i)] .+= Ub[i][1]
        Ubvec.v[getrange(bodies,i)] .+= Ub[i][2]
    end
end

function _computef̃limit(SP::SuctionParameter, plate::Plate, Γ₀)
    f̃ = -SP*2π*plate.len/Γ₀
    return f̃
end

function _computef̃limits(SP::SuctionParameter, args...)
    return _computef̃limits(SuctionParameterRange(-SP,SP), args...)
end

function _computef̃limits(SP::SuctionParameterRange, args...)
    f̃min = _computef̃limit(SP.σmax, args...)
    f̃max = _computef̃limit(SP.σmin, args...)
    return f̃Limits(f̃min,f̃max)
end
