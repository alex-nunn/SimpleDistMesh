import Base.show
using LinearAlgebra
using ForwardDiff
using DiffResults
using StaticArrays

"""
    SignedDistFunc

Abstract type for all signed distance functions
"""
abstract type SignedDistFunc <: Function end


# ------------------------------------------------------------------------------
# Combining signed distance functions
# ------------------------------------------------------------------------------
"""
`union(f, g)` returns the union of the regions for signed distance functions 
"""
Base.union(f1::Function, f2::Function) = pt -> min(f1(pt), f2(pt))

"""
`setdiff(f, g)` returns the region of `f` with `g` removed for signed distance
functions
"""
Base.setdiff(f1::Function, f2::Function) = pt -> max(f1(pt), -f2(pt))

"""
`intersect(f, g)` returns the intersection of the `f` and `g` regions for signed
distance functions
"""
Base.intersect(f1::Function, f2::Function) = pt -> max(f1(pt), f2(pt))


# ------------------------------------------------------------------------------
# Circle
# ------------------------------------------------------------------------------
"""
    Circle(p0, r)

Signed distance function for a circle with center `p0` and radius `r`
"""
struct Circle <: SignedDistFunc
    p0
    r
end
(c::Circle)(pt) = norm(pt - c.p0) - c.r
Base.show(io::IO, ::MIME"text/plain", c::Circle) = print(io, "Circle($(c.p0), $(c.r))")


# ------------------------------------------------------------------------------
# Rectangle
# ------------------------------------------------------------------------------
"""
    Rect(x1, x2, y1, y2)

Signed distance function for a rectangle `[x1, x2] × [y1, y2]`
"""
struct Rect <: SignedDistFunc
    x1
    x2
    y1
    y2
end
(r::Rect)(pt) = -min(-r.y1+pt[2], r.y2-pt[2], -r.x1 + pt[1], r.x2 - pt[1])
Base.show(io::IO, ::MIME"text/plain", r::Rect) = print(io, "Rect([$(r.x1), $(r.x2)]×[$(r.y1), $(r.y2)])")


# ------------------------------------------------------------------------------
# Implicit region
# ------------------------------------------------------------------------------
"""
    ImplicitRegion(f)

Signed distance function for the implicit region defined by `f(x) = 0`
"""
struct ImplicitRegion <: SignedDistFunc
    f::Function
end

"""
`ImplicitRegion(f, c)` returns a signed distance function for the implicit
region defined by `f(x) = c`
"""
ImplicitRegion(f, c) = ImplicitRegion(x -> f(x) - c)
(region::ImplicitRegion)(x0) = begin
    x = implicit_project(region.f, x0)
    return sign(region.f(x0)) * norm(x - x0)
end

"""
    implicit_project(f, x0; max_iters, α)

Return the projection of `x0` onto the level set `f(x) = 0` such that the
displacement from `x0` is parallel to the gradient, `(x - x0) ⋅ ∇f(x) = 0`
"""
function implicit_project(f, x0; max_iters=1000, α=1.0)
    # Allocate memory
    L = zeros(MVector{2})
    Δx = zeros(MVector{2})
    J = zeros(MMatrix{2, 2})
    H = zeros(MMatrix{2, 2})
    x = zeros(2)
    dr = DiffResults.HessianResult(x0)
    hconf = ForwardDiff.HessianConfig(f, dr, x0)

    x .= x0
    for i ∈ 1:max_iters
        ForwardDiff.hessian!(dr, f, x, hconf)
        ∇fx, ∇fy = DiffResults.gradient(dr)

        L[1] = DiffResults.value(dr)
        L[2] = (x[1] - x0[1]) * ∇fy - (x[2] - x0[2]) * ∇fx

        H .= DiffResults.hessian(dr)
        J[1, 1] = ∇fx
        J[1, 2] = ∇fy
        J[2, 1] = ∇fy + (x[1] - x0[1]) * H[1, 2] - (x[2] - x0[2]) * H[1, 1]
        J[2, 2] = -∇fx - (x[2] - x0[2]) * H[2, 1]  + (x[1] - x0[1]) * H[2, 2]

        Δx .= inv(J) * L
        @. x -= α * Δx
        if (Δx[1]^2 + Δx[2]^2) < 1e-8
            break
        end
    end
    return x
end


# ------------------------------------------------------------------------------
# Polygon
# ------------------------------------------------------------------------------
"""
    Polygon(nodes)
    
Signed distance function for a polygon with 2-dimensional nodes as columns of a
matrix

# Examples
Create signed distance function for a triangle
```
julia> triangle = Polygon([0 1 0; 0 0 1])
Polygon([0, 0], [1, 0], [0, 1])
```
"""
struct Polygon <: SignedDistFunc
    nodes::Matrix
end
(poly::Polygon)(x::AbstractVector) = (x ∈ poly ? -1 : 1) * distance(poly, x)
Base.show(io::IO, ::MIME"text/plain", p::Polygon) = begin
    print(io, "Polygon($(join(eachcol(p.nodes), ", ")))")
end

"""
    distance(poly, x)

Return the minimal 2-norm distance from `x` to the edges of the polygon `poly`
"""
function distance(poly::Polygon, x::AbstractVector)
    d = Inf
    a = poly.nodes[:, end]
    for b ∈ eachcol(poly.nodes)
        t = (b - a)⋅(x - a) / ((b[1] - a[1])^2 + (b[2] - a[2])^2)

        if t <= 0
            dp = norm(x - a)
        elseif t >= 1
            dp = norm(x - b)
        else
            dp = norm(a - x + (b - a) * t)
        end

        if dp < d
            d = dp
        end
        a = b
    end
    return d
end

"""
    winding_number(poly, x)

Return the winding number of polygon `poly` about point `x`
"""
function winding_number(poly::Polygon, x::AbstractVector)
    Γ = 0  # winding number about point x

    a = poly.nodes[:, end]
    for b ∈ eachcol(poly.nodes)
        if a[2] <= x[2]
            if x[2] < b[2] && triangle_area(a, b, x) > 0
                Γ += 1
            end
        elseif b[2] <= x[2] && triangle_area(a, b, x) < 0
            Γ -= 1
        end
        a = b
    end
    return Γ
end

"""
    in(x::Vector, poly::Polygon)

Return true if point `x` is inside polygon `poly`
"""
Base.in(x::AbstractVector, poly::Polygon) = isodd(winding_number(poly, x))

"""
    triangle_area(a, b, c)

Return the signed area of the 2D triangle `abc`
"""
@inline function triangle_area(a, b, c)
    return 0.5((a[1] - c[1]) * (b[2] - c[2]) - (a[2] - c[2]) * (b[1] - c[1]))
end


# ------------------------------------------------------------------------------
# Rotations
# ------------------------------------------------------------------------------
"""
    Rotation(ϕ, p0)

Rotation about point `p0` counter-clockwise by angle `ϕ` (in radians)

# Examples
Rotate point about the origin by `π/2`
```
julia> Rotation(π/2)([1.0, 0.0])
2-element Vector{Float64}:
 6.123233995736766e-17
 1.0
```
"""
struct Rotation <: Function
    A::Matrix
    p0
    Rotation(ϕ, p0=[0.0, 0.0]) = new([cos(ϕ) -sin(ϕ); sin(ϕ) cos(ϕ)], p0)
end
(rot::Rotation)(pt::Vector) = rot.A * (pt - rot.p0) + rot.p0 