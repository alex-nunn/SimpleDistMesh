# SimpleDistMesh.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://alex-nunn.github.io/SimpleDistMesh/)

A Julia implementation of [DistMesh](http://persson.berkeley.edu/distmesh/), a simple signed distance function meshing algorithm developed by [Per-Olof Persson](http://persson.berkeley.edu/).

One motivation for the DistMesh algorithm is to provide a simple meshing algorithm, easily understood and modified by newcomers. Given the importance of mesh generation across a broad range of discplines, _"the ability to understand and adapt mesh generation code is too valuable an option to lose."_ Thus, whilst the algorithm is still capable of producing high quality meshes, performance is sacrificed for clarity.[1]

Some minor changes have been made to the MATLAB implementation to express the ideas more idiomatically in Julia.

## References
[1] [Persson, P. O., & Strang, G. (2004). A simple mesh generator in MATLAB. SIAM review, 46(2), 329-345.](https://doi.org/10.1137/S0036144503429121)

[2] [Persson, P. O. (2005). Mesh generation for implicit geometries (Doctoral dissertation, Massachusetts Institute of Technology).](http://dspace.mit.edu/handle/1721.1/27866)