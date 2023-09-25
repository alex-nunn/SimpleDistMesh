# Docs make.jl

push!(LOAD_PATH, "../src/")
using SimpleDistMesh
using Documenter

makedocs(
    sitename="SimpleDistMesh.jl",
    modules=[SimpleDistMesh],
    pages=[
        "Home"=>"index.md",
        "Examples"=>"examples.md",
        "Functions"=>[
            "Meshing"=>"meshing.md",
            "Signed distance functions"=>"signed_distance_functions.md"
        ]
    ],
)

deploydocs(;
    repo="github.com/alex-nunn/SimpleDistMesh.jl"
)