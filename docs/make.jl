# Docs make.jl

push!(LOAD_PATH, "../src/")
using SimpleDistMesh
using Documenter

makedocs(
    sitename="SimpleDistMesh.jl",
    modules=[SimpleDistMesh],
    pages=[
        "Home"=>"index.md"
    ]
)

#deploydocs(;
#    repo="github.com/alex-nunn/SimpleDistMesh.jl"
#)