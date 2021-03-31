# GridPotentialFlow.jl

*A set of tools to solve potential flows past bodies on a Cartesian grid.*

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4549940.svg)](https://doi.org/10.5281/zenodo.4549940)

The objective of this package is to allow easy setup and fast simulation of potential
flows. The package provides tools for
- constructing grids, body shapes, and point vortices,
- specifying the relevant parameters and setting their values,
- specifying the edges where shedding occurs and setting their suction parameter,
- solving the problem.

The underlying grids are uniform and Cartesian, making use of the [CartesianGrids](https://github.com/JuliaIBPM/CartesianGrids.jl) package. This package allows the use of the lattice Green's function (LGF) for inverting the Poisson equation. The presence of bodies is accounted for using the immersed boundary projection method, originally developed for viscous flow by Taira and Colonius [^1]. The potential flow system with the no-penetration condition, any edge conditions, and  Kelvin's circulation theorem is implemented with the [ConstrainedSystems](https://github.com/JuliaIBPM/ConstrainedSystems.jl) package. Tools for creating bodies are based on the [RigidBodyTools](https://github.com/JuliaIBPM/RigidBodyTools.jl) package. The vortex dynamics are computed using the vortex-in-cell method of Christiansen [^2]. For more details, please refer to [Beckers, D. and Eldredge, J. D. (2021) "Planar potential flow on Cartesian grids," [arXiv:2102.11910]](https://arxiv.org/abs/2102.11910).

![vortexshedding](https://user-images.githubusercontent.com/26737762/113199963-9b77ee80-921c-11eb-8448-70a32e50660f.gif)

The plots in this documentation are generated using [Plots.jl](http://docs.juliaplots.org/latest/).
You might want to install that, too, to follow the examples.

## References

[^1]: Taira, K. and Colonius, T. (2007) "The immersed boundary method: a projection approach," *J. Comput. Phys.*, 225, 2118--2137.

[^2]: Christiansen, J. (1973) "Numerical simulation of hydrodynamics by the method of point vortices," *J. Comput. Phys.*, 13, 363--379.
