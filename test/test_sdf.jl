using Test
using SimpleDistMesh.sdf

@testset "Polygon" begin
    psqr = Polygon([
        0.0 1.0 1.0 0.0
        0.0 0.0 1.0 1.0
    ])

    @testset "triangle_area" begin
        @test sdf.triangle_area([0.0, 0.0], [1.0, 0.0], [0.0, 1.0]) ≈ 0.5
    end

    @testset "winding_number" begin
        @test sdf.winding_number(psqr, [0.5, 0.5]) == 1
        @test sdf.winding_number(psqr, [-0.5, 0.5]) == 0
        @test sdf.winding_number(psqr, [1.5, 0.5]) == 0
    end

    @testset "Base.in" begin
        @test [0.5, 0.5] ∈ psqr
        @test [-0.5, 0.5] ∉ psqr
        @test [1.5, 0.5] ∉ psqr
    end

    @testset "distance" begin
        @test sdf.distance(psqr, [-0.5, 0.1]) ≈ 0.5
    end
end;