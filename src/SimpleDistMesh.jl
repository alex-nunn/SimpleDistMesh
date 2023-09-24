"""
    SimpleDistMesh

Mesh generation using signed distance functions

# Examples
Create a mesh for a circular disk
```
julia> mesh = Mesh(
           Circle([0, 0], 1),  # d -- signed distance function
           pt -> 1,            # h -- element size function
           0.05,               # h0 -- initial element size
           [-1 -1; 1, 1]       # bounds -- [x_min y_min; x_max y_max]
       )
```

Mesh for a square `[-1, 1]×[-1, 1]` with a circular hole
```
julia> d = setdiff(
           Rect(-1, 1, -1, 1)),
           Circle([0.0, 0.0], 1/2)
       );
julia> mesh = Mesh(
           d,                  # d -- signed distance function
           pt -> 1,            # h -- element size function
           0.05,               # h0 -- initial element size
           [-1 -1; 1, 1]       # bounds -- [x_min y_min; x_max y_max]
       );
```
"""
module SimpleDistMesh

include("distmesh.jl")
include("sdf.jl")

end # module SimpleDistMesh
