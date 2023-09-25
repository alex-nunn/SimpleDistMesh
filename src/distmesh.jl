"""
    Mesh

Mesh with node positions and triangulation `ts`

The columns of the `nodes` matrix provide  the positions of points in the mesh, 
while the triangles (or tetrahedra) of the mesh defined by the columns of `ts`.

# Examples
A mesh in 2D comprised of the single triangle (0, 0), (1, 0), (1, 1)

julia> mesh = Mesh{Float64, 2}(
            [ 
                0.0 1.0 0.0
                0.0 0.0 1.0
            ],
            [1; 2; 3]
        )
"""
struct Mesh
    nodes::Matrix  # positions of nodes
    ts::Matrix     # triangulation
end


"""
    Mesh(d, h, h0, bounds, fixed_nodes; ...)

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
- `Fscale` internal pressure of relaxation method
- `fixed_nodes::Matrix` defines fixed points of the mesh (if any)
- `Δt` time-step of Euler method for point relaxation
- `geps` tolerance in geometry calculations
- `dptol` stopping criterion for relaxation method. Terminate if all node 
    advections move less than `dptol`.
- `ttol` for computational efficiency, triangulations are not computed at every
    step of the mesh relaxation. Recompute the mesh triangulation if a point has
    moved by at least `ttol` since last triangulation.
- `max_iters` maximum number of iterations

# Outputs
Returns a `Mesh`
"""
function Mesh(
        d::Function, 
        h::Function, 
        h0::Real,
        bounds;
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
    bars, ts = find_bars(d, geps, nodes) # Compute bar connections between nodes

    for i ∈ 1:max_iters
        # Retriangulation by the Delaunay algorithm
        # (Recompute bars if nodes exceed movement tolerance ttol)
        Δ = maximum(norm.(eachcol(nodes - prev_nodes))) / h0
        if Δ > ttol
            prev_nodes[:] .= nodes[:]
            bars, ts = find_bars(d, geps, nodes)
        end

        # Update node positions 
        δinterior = relax_mesh!(nodes, bars; d, h, n_fixed_nodes, fscale, Δt)

        # Termination criterion
        if δinterior < dptol * h0
            println("Successfully converged after $i iterations.")
            return Mesh(nodes, ts)
        end
    end

    println("Maximum number of iterations reached. Terminating.")

    return Mesh(nodes, ts)
end


"""
    relax_mesh!(nodes, ts; ...)

Perform one relaxation iteration on the `nodes` and mesh triangulation `ts and
return the maximum distance moved by an iterior node.
"""
function relax_mesh!(
        nodes::Matrix{T}, bars; 
        d, h, n_fixed_nodes, fscale, Δt
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
            nodes[:, i] -= d0 * ForwardDiff.gradient(d, node)
        else
            δ = Δt * norm(node_forces[:, i])
            δinterior = max(δinterior, δ)
        end
    end
    return δinterior
end


"""
    grid_uniform_triangles(bounds, h0)

Return nodes of a uniform grid of equilaterial triangles with side-length `h0`. 
Nodes are returned as a vector of vectors.
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

Return subject of nodes inside the domain which pass the rejection method
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

Return the bars which connect the `nodes` of the mesh using Delaunay 
triangulation and geometrical tolerance `geps`. The bars are defined by the 
columns the returned matrix.
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


@recipe function recipe_mesh(mesh::Mesh)
    aspect_ratio --> :equal

    @series begin
        label --> false
        seriestype := :scatter
        (mesh.nodes[1, :], mesh.nodes[2, :])
    end

    for tri ∈ eachcol(mesh.ts)
        tri = [tri..., tri[1]]
        @series begin
            label --> false
            seriestype := :path
            fillcolor --> false
            linewidth --> 1
            linecolor --> palette(:tab10)[begin]

            mesh.nodes[1, tri], mesh.nodes[2, tri]
        end
    end
end
