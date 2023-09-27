"""
    SimpleDistMesh

Mesh generation using signed distance functions
"""
module SimpleDistMesh

using RecipesBase
using DelaunayTriangulation
using ForwardDiff
using LinearAlgebra

export Mesh
export Circle, Rect, Rotation, Polygon, ImplicitRegion

include("distmesh.jl")
include("sdf.jl")

end # module SimpleDistMesh
