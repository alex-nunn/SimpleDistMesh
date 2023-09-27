# ------------------------------------------------------------------------------
# Build assets for examples documentation
# ------------------------------------------------------------------------------
using SimpleDistMesh
using Plots
using Plots.Measures
using LinearAlgebra

# Setup export 
assets_dir = joinpath(Base.source_dir(), "src", "assets")
fig_ouput = (p, name) -> savefig(p, joinpath(assets_dir, name))
mkpath(assets_dir)

# Plot defaults
pkws = (
    markersize=2,
    margin=0mm,
    frame=:none,
    dpi=200
)

group_plot = meshes -> plot(
    plot.(meshes)...;
    layout=(1, length(meshes)), 
    size=(length(meshes) * 330, 330), 
    pkws...
)


# ------------------------------------------------------------------------------
# 1. Disks of increasing discretisation
# ------------------------------------------------------------------------------
d = Circle([0.0, 0.0], 1)
h = x -> 1.0
bounds = [-1 -1; 1 1]

h0s = [0.4, 0.2, 0.1]
meshes = [Mesh(d, h, h0, bounds) for h0 ∈ h0s]

p = group_plot(meshes)
fig_ouput(p, "01_disks.png")

# ------------------------------------------------------------------------------
# 2. Combining regions
# ------------------------------------------------------------------------------
# 2.1 Circular annulus
c_outer = Circle([0.0, 0.0], 1)
c_inner = Circle([0.0, 0.0], 0.3)
d = setdiff(c_outer, c_inner)

h = x -> 1.0
h0 = 0.1
bounds = [-1 -1; 1 1]

mesh1 = Mesh(d, h, h0, bounds)

# 2.2 Square with hole
square = Rect(-1, 1, -1, 1)
hole = Circle([0.0, 0.0], 0.3)
d = setdiff(square, hole)

h = x -> 1.0
h0 = 0.15
bounds = [-1 -1; 1 1]
fixed_nodes = [-1 -1 1 1; -1 1 -1 1]

mesh2 = Mesh(d, h, h0, bounds; fixed_nodes)

# 2.3 Square with hole (adaptive sizing)
square = Rect(-1, 1, -1, 1)
hole = Circle([0.0, 0.0], 0.4)
d = setdiff(square, hole)

h = x -> min(0.2 + norm(x) - 0.4, 1)
h0 = 0.05
bounds = [-1 -1; 1 1]
fixed_nodes = [-1 -1 1 1; -1 1 -1 1]

mesh3 = Mesh(d, h, h0, bounds; fixed_nodes)

meshes = [mesh1, mesh2, mesh3]
p = group_plot(meshes)
fig_ouput(p, "02_combining_regions.png")

# ------------------------------------------------------------------------------
# 3. Complex geometries
# ------------------------------------------------------------------------------
# 3.1 Hexagon
hexagon = (;l, ϕ0) -> Polygon(
    stack((l * cos(ϕ), l * sin(ϕ)) for ϕ ∈ ϕ0 .+ range(0, 2π, 7)[begin:end-1])
)

outer_hexagon = hexagon(l=1.0, ϕ0=0)
inner_hexagon = hexagon(l=0.5, ϕ0=π/6)

d = setdiff(outer_hexagon, inner_hexagon)
h = x -> 1
h0 = 0.1
bounds = transpose(stack(extrema(outer_hexagon.nodes; dims=2);dims=1))
fixed_nodes = inner_hexagon.nodes

mesh1 = Mesh(d, h, h0, bounds; fixed_nodes)

# 3.2 Geometric adaptivity
upper_half_plane = x -> -x[2]
d1 = Circle([0, 0], 1)
d2 = Circle([-0.4, 0.0], 0.55)
d = setdiff(d1, d2) ∩ upper_half_plane

h1 = x -> 0.15 - 0.2d1(x)
h2 = x -> 0.06 + 0.2d2(x)
h3 = x -> (d2(x) - d1(x)) / 3

h = x -> min(h1(x), h2(x), h3(x))
h0 = 0.05 / 3
bounds = [-1 0; 1 1]
fixed_nodes = [
    -1.0 -0.95 0.15 1.0
     0.0  0.00 0.00 0.0
]

mesh2 = Mesh(d, h, h0, bounds; fixed_nodes)

meshes = [mesh1, mesh2]
p = group_plot(meshes)
fig_ouput(p, "03_complex_geometries.png")


# ------------------------------------------------------------------------------
# 4. Implicit regions
# ------------------------------------------------------------------------------
# 4.1 Super ellipses
d = setdiff(
    ImplicitRegion(x -> norm(x, 4) - 1),
    ImplicitRegion(x -> norm(x, 4) - 0.5)
)

h = x -> 1
h0 = 0.1
bounds = [-1 -1; 1 1]


mesh1 = Mesh(d, h, h0, bounds)

# 4.2 Weird region
d = ImplicitRegion(x -> x[2] - cos(x[1])) ∩ ImplicitRegion(x -> -x[2] + 5((2 * x[1] / (5π))^4 - 1))
h = x -> 1
h0 = 0.6
bounds = [-5π/2 -5; 5π/2 1]
fixed_nodes = [-5π/2 5π/2; 0 0]

mesh2 = Mesh(d, h, h0, bounds; fixed_nodes)

meshes = [mesh1, mesh2]
p = group_plot(meshes)
fig_ouput(p, "04_complex_geometries.png")

# ------------------------------------------------------------------------------
# 5. All together now
# ------------------------------------------------------------------------------
r = norm([-4, 0] - [1.0, 0.5])
ϕ = atan(1.5, 5)

d1 = Polygon([
    -4.0  1.0  1.0  0.0  1.0  1.0
     0.0  1.5  0.5  0.0 -0.5 -1.5
]) ∩ Circle([-4, 0], r)
d2 = Circle([-1.0, 0.0], 0.5)


d = setdiff(d1, d2)
h1 = x -> clamp(0.7norm(x), 0.1, 0.8)
h2 = x -> 0.2 + d2(x)
h = x -> min(h1(x), h2(x))

h0 = 0.03
bounds = [-4 -1.5; 1 1.5]

fixed_nodes = stack([
    [-4.0, 0.0],
    [1.0, 0.5],
    [1.0, -0.5],
    [0.0, 0.0],
    [r * cos(ϕ) - 4, r * sin(ϕ)],
    [r * cos(ϕ) - 4, -r * sin(ϕ)]
])

mesh = Mesh(d, h, h0, bounds; fixed_nodes)
p = plot(
    mesh;
    size=(600, 400),
    markersize=0,
    margin=0mm,
    frame=:none,
    dpi=200
)
fig_ouput(p, "05_all_together_now.png")
