function computeregularizationmatrix(g::PhysicalGrid, X::VectorData{N}, f::ScalarData{N}, s::Nodes; ddftype=CartesianGrids.Yang3) where {N}

    regop = Regularize(X, cellsize(g), I0=origin(g), ddftype = ddftype, issymmetric=true)
    Rmat,_ = RegularizationMatrix(regop, f, s)
    Emat = InterpolationMatrix(regop, s, f)

    return Rmat, Emat
end

function computecrossproductsurfaceintegral(body::Body, z)
    rx,ry = collect(body)
    @assert length(rx) == length(z)
    return [ry'*z, -rx'*z]
end

function _computebodypointsvelocity!(Ubvec,Ub,bodies)

    # Assure that Ub is an array
    if Ub isa Tuple
        Ub = [Ub]
    end

    @assert length(Ub) == length(bodies)

    Ubvec.u .= 0
    Ubvec.v .= 0

    for i in 1:length(bodies)
        Ubvec.u[getrange(bodies,i)] .+= Ub[i][1]
        Ubvec.v[getrange(bodies,i)] .+= Ub[i][2]
    end
end
