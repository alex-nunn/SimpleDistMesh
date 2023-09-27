```@index
Pages=["mesh.md"]
```

```@meta
CurrentModule = SimpleDistMesh
```

# Mesh
Methods for creating 2-dimensional meshes using the `DistMesh` method are contained in the module,
```@docs
SimpleDistMesh
```

---

Meshes are stored as a set of node positions and triangles in column-major format in `Mesh`. A 2-dimensional mesh can be generated by defining a signed distance function $d:\mathbb{R}^2 \rightarrow  \mathbb{R}$ and desired truss element length $h:\mathbb{R}^2 \rightarrow  \mathbb{R}$, and passing to the `Mesh` constructor.
```@docs
Mesh
```

---
### Utilities

The following functions perform various steps of the `DistMesh` mesh relaxation algorithm.

```@docs
relax_mesh!
grid_uniform_triangles
rejection_method
find_bars
```