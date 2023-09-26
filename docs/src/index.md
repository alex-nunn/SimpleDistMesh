# SimpleDistMesh.jl

A Julia implementation of [DistMesh](http://persson.berkeley.edu/distmesh/), a simple signed distance function meshing algorithm developed by [Per-Olof Persson](http://persson.berkeley.edu/).

 Given the importance of mesh generation across a broad range of discplines, _"the ability to understand and adapt mesh generation code is too valuable an option to lose."_ [1] Thus, a motivation for the DistMesh algorithm is to provide a simple meshing algorithm capable of producing high quality meshes, easily understood and modified by newcomers. 
 
In keeping with this goal our implementation sacrifices performance in favor of readability. Also, some minor changes have been made to the MATLAB implementation to express the ideas more idiomatically in Julia.

## References
[1] [Persson, P. O., & Strang, G. (2004). A simple mesh generator in MATLAB. SIAM review, 46(2), 329-345.](https://doi.org/10.1137/S0036144503429121)

[2] [Persson, P. O. (2005). Mesh generation for implicit geometries (Doctoral dissertation, Massachusetts Institute of Technology).](http://dspace.mit.edu/handle/1721.1/27866)