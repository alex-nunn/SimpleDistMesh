"""
    SimpleDistMesh

Mesh generation using signed distance functions
"""
module SimpleDistMesh

export Mesh
export Circle, Rect, Rotation, Polygon, ImplicitRegion

include("distmesh.jl")
include("sdf.jl")

end # module SimpleDistMesh
