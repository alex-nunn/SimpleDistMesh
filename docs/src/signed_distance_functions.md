```@index
Pages=["signed_distance_functions.md"]
```
```@meta
CurrentModule = SimpleDistMesh
```
# Signed distance functions
A signed distance function, $d(x)$, can be any scalar function in 2-dimensions for which $d(x) = 0$ defines a meaningful region. For convenience, we have provided several simple types for creating common shapes.
```@docs
Circle
Rect
Polygon
```

```@example
using SimpleDistMesh
c = Circle([0.0, 0.0], 1.0)
```

Signed distance functions can also be composed with coordinate transforms to define a range of different geoemtries. For convenience, we also provide a type for defining rotations.
```@docs
Rotation
```

All of these are subtypes of the abstract type,
```@docs
SignedDistFunc
```
---

The regions represented by signed distance can also be combined using common set operations of unions, intersections and set differences.

```@docs
union
setdiff
intersect
```

---
### Utilities
To create signed distance functions for polygonal geometries we define several useful functions.

```@docs
distance
winding_number
in
triangle_area
```