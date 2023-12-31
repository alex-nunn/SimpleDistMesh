using FiniteDiff
using RecipesBase
using DelaunayTriangulation
using ForwardDiff
using LinearAlgebra

"""
    Mesh

Mesh containing `nodes` positions and triangulation

# Properties
- `nodes::Matrix` (2, n_nodes) 2-dimensional position of nodes
- `triangulation::Matrix{Int}` (3, n_edges) node indices defining the triangles
    of the mesh
"""
struct Mesh
    nodes::Matrix             # positions of nodes
    triangulation::Matrix     # connections between nodes
end


"""
    Mesh(d, h, h0, bounds; ...)

Construct a mesh using the signed distance function `d` and desired edge length
function `h`. 
    
The initial edge length is given by `h0` and the mesh is constructed within the 
bounding box given by `bounds`. Fixed nodes of the mesh can be provided as 
columns of `fixed_nodes`.

# Arguments
- `d::Function` signed distance function (== 0 on the boundary)
- `h::Function` preferred edge length function
- `h0` starting edge length
- `bounds::Matrix` defines bounding box for mesh with columns defining the 
    minimum and maximum for each axis. For example,
    `bounds = [x_min y_min; x_max y_max]`

# Keyword Arguments
- `gradient` function with signature `gradient(d, x)` returning the gradient at
        x (default: `FiniteDiff.finite_difference_gradient`)
- `fixed_nodes::Matrix` defines fixed points of the mesh (if any)
- `Fscale` internal pressure of relaxation method (default: 1.2)
- `Δt` time-step of Euler method for point relaxation (default: 0.2)
- `geps` tolerance in geometry calculations (default: 0.001 h0)
- `dptol` stopping criterion for relaxation method. Terminate if all nodes 
    move less than `dptol` in a single iteration. (default: 0.001)
- `ttol` for computational efficiency, triangulations are not computed at every
    step of the mesh relaxation. Recompute the mesh triangulation if a point has
    moved by at least `ttol` since last triangulation. (default: 0.1)
- `max_iters` maximum number of iterations (default: 1000)

# Outputs
Returns a `Mesh`
"""
function Mesh(
        d::Function, 
        h::Function, 
        h0::Real,
        bounds;
        gradient=FiniteDiff.finite_difference_gradient,
        fixed_nodes=Matrix{Float64}(undef, 2, 0),
        fscale::Real=1.2,
        Δt::Real=0.2,
        geps::Real=0.001 * h0,
        dptol::Real=0.001,
        ttol::Real=0.1,
        max_iters=10^3
        )

    # Create a uniform grid of equilateral triangles
    nodes = grid_uniform_triangles(bounds, h0)

    # Remove points outside of the domain and apply rejection method
    nodes = rejection_method(d, h, geps, nodes)
    nodes = unique(hcat(fixed_nodes, nodes); dims=2)  # add in fixed nodes
    n_fixed_nodes = size(fixed_nodes, 2)
    
    # Relaxation procedure
    prev_nodes = copy(nodes)

    # Compute bar connections between nodes
    bars, triangulation = find_bars(d, geps, nodes) 

    for i ∈ 1:max_iters
        # Retriangulation by the Delaunay algorithm
        # (Recompute bars if nodes exceed movement tolerance ttol)
        Δ = maximum(norm.(eachcol(nodes - prev_nodes))) / h0
        if Δ > ttol
            prev_nodes[:] .= nodes[:]
            bars, triangulation = find_bars(d, geps, nodes)
        end

        # Update node positions 
        δinterior = relax_mesh!(
            nodes, bars; d, h, n_fixed_nodes, gradient, fscale, Δt
        )

        # Termination criterion
        if δinterior < dptol * h0
            println("Successfully converged after $i iterations.")
            return Mesh(nodes, triangulation)
        end
    end

    println("Maximum number of iterations reached. Terminating.")

    return Mesh(nodes, triangulation)
end


"""
    relax_mesh!(nodes, bars; ...)

Move the nodes of the mesh according the forces between nodes, and return the 
maximum distance moved by an iterior node.
"""
function relax_mesh!(
        nodes::Matrix{T}, bars; 
        d, h, n_fixed_nodes, gradient, fscale, Δt
        ) where T <: Real
    
    # 6. Move mesh points bar on bar lengths Ls and forces `node_forces`
    n_bars = size(bars, 2)
    Ls = Vector{T}(undef, n_bars)       # bar lengths
    hbars = Vector{T}(undef, n_bars)    # size function at bar midpoints

    # Loop through bars
    for i ∈ axes(bars, 2)
        n1, n2 = bars[:, i]
        p1 = nodes[:, n1]
        p2 = nodes[:, n2]

        Ls[i] = norm(p2 - p1)
        hbars[i] = h(0.5(p1 + p2))
    end

    # Desired lengths
    L0s = hbars * fscale * sqrt(sum(Ls .^ 2) / sum(hbars .^ 2))

    # Empty matrix for force on each node
    node_forces = zeros(T, size(nodes)...)

    # Loop through bars
    for i ∈ axes(bars, 2)
        n1, n2 = bars[:, i]
        p1 = nodes[:, n1]
        p2 = nodes[:, n2]
        L = Ls[i]
        L0 = L0s[i]

        force = max(L0 - L, 0) * (p1 - p2) / L
        node_forces[:, n1] .+= force
        node_forces[:, n2] .-= force
    end
    
    node_forces[:, 1:n_fixed_nodes] .= 0.0  # set movement of fixed nodes to zero

    # Update node positions using Euler's method
    @. nodes += Δt * node_forces

    # Project outside points back onto the boundary
    δinterior = 0.0
    for i ∈ axes(nodes, 2)
        node = nodes[:, i]
        d0 = d(node)

        # Is node outside domain?
        if d0 > 0
            nodes[:, i] -= d0 * gradient(d, node)
        else
            δ = Δt * norm(node_forces[:, i])
            δinterior = max(δinterior, δ)
        end
    end
    return δinterior
end


"""
    grid_uniform_triangles(bounds, h0)

Return the nodes for a grid of equilaterial triangles with side-length `h0`. 
"""
function grid_uniform_triangles(bounds, h0)
    xs = range(bounds[:, 1]...; step=h0)
    ys = range(bounds[:, 2]...; step=h0 * sqrt(3) / 2)
    nodes = Matrix{Float64}(undef, 2, length(xs) * length(ys))

    nidx = 1
    for (i, x) ∈ enumerate(xs), y ∈ ys
        nodes[1, nidx] = x
        nodes[2, nidx] = iseven(i) ? y + 0.5h0 : y
        nidx += 1
    end
    return nodes
end


"""
    rejection_method(d, h, nodes)

Return interior nodes which pass the rejection method
"""
function rejection_method(d, h, geps, nodes)
    idxs = []
    r0s = []
    for (i, node) ∈ enumerate(eachcol(nodes))
        if d(node) < geps
            push!(idxs, i)  # record index of node inside domain
            push!(r0s, 1 / h(node)^2)  # compute rejection factor
        end
    end
    r0s /= maximum(r0s)  # divide by maximum to obtain probability

    idxs = idxs[rand(length(idxs)) .< r0s]  # reject nodes
    nodes = nodes[:, idxs]
    return nodes
end


"""
    find_bars(d, geps, nodes)

Using a Delaunay triangulation return the bars connecting the `nodes` of the 
mesh and the indices of the mesh triangles

# Output
`(bars, triangulation)`

- `bars::Matrix{Int}` (2, n_edges) indices of nodes defining bar connections
- `triangulation::Matrix{Int}` (3, n_triangles) indices of nodes defining the
    triangles of the mesh
"""
function find_bars(d, geps, nodes)
    tri = triangulate(nodes)  # compute Delaunay triangulation
    bar_set = Set{NTuple{2, Int}}()
    tri_set = Vector{NTuple{3, Int}}()

    for t ∈ each_solid_triangle(tri)
        t1, t2, t3 = sort(t)
        p1 = nodes[:, t1]
        p2 = nodes[:, t2]
        p3 = nodes[:, t3]
        
        centroid = (p1 + p2 + p3) / 3
        # Is centroid within domain interior?
        if d(centroid) < -geps
            push!(tri_set, t)
            push!(bar_set, (t1, t2), (t1, t3), (t2, t3))
        end
    end
    return stack(bar_set), stack(tri_set)
end

# Recpie for plotting mesh objects
@recipe function recipe_mesh(mesh::Mesh)
    aspect_ratio --> :equal

    for tri ∈ eachcol(mesh.triangulation)
        tri = [tri..., tri[1]]
        @series begin
            label --> false
            seriestype := :path
            fillcolor --> false
            linewidth --> 1
            linecolor --> "#1F77B4"

            mesh.nodes[1, tri], mesh.nodes[2, tri]
        end
    end

    @series begin
        label --> false
        seriestype := :scatter
        markercolor --> "#1F77B4"

        (mesh.nodes[1, :], mesh.nodes[2, :])
    end
end
