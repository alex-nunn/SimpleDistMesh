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
export make_mesh, relax_mesh!, grid_uniform_triangles, rejection_method

include("distmesh.jl")
include("sdf.jl")

end # module SimpleDistMesh
