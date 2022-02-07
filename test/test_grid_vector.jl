using WTP
using Test
using LinearAlgebra
using Profile
using StaticArrays

@testset "Grid Vector" begin

    # Create a grid vector
    reciprocal_basis = CARTESIAN_BASIS
    sizes = (4, 4, 4)
    lattice = make_grid(ReciprocalLattice3D, reciprocal_basis, size_to_domain(sizes))
    g = GridVector{ReciprocalLattice3D}(lattice, [0, 1, 2], true)

    # Get and set grid.
    @test grid(g) == lattice
    # THe size of the grid vector should be that of the grid.
    @test size(grid(g)) == sizes
    new_g = set_grid(g, expand(lattice, [1, 2, 3]))
    @test grid(new_g) != grid(g)

    # Indexing a grid should yield a grid vector.
    @test lattice[0, 1, 2] == g

    # Wrapping.
    h = lattice[0, 0, 3]
    @test wrapped(h) == [0, 0, -1]
    @test wrapped(lattice[0, 0, -3]) == [0, 0, 1]
    @test wrapped(lattice[3, -3, -3]) == [-1, 1, 1]

    # Overflow
    @test overflow(h) == [0, 0, 1]
    @test overflow(lattice[0, 0, -3]) == [0, 0, -1]
    @test overflow(lattice[3, -3, -3]) == [1, -1, -1]
    # Wrapping should compute with negation.
    @test has_overflow(h)

    # Get the image within the grid domain.
    h2 = reset_overflow(h)
    @test !has_overflow(h2)

    # Get back the image out of the grid domain.
    @test h == add_overflow(h2, overflow(h))
    @test 9 == braket(dagger(h), h)
    # Snapping after taking the coordinates should recover the vector
    @test snap(lattice, [0, 1, 2]) == g

    # Basic vector algebra.
    @test add(lattice[0, 0, 1], lattice[0, 1, 0]) == lattice[0, 1, 1]
    # Mismatching grids.
    @test_broken add(lattice[0, 0, 1], make_grid_vector(expand(lattice), SVector(0, 1, 0)))
    @test negate(lattice[0, 0, 1]) == lattice[0, 0, -1]
    @test mul(2, lattice[1, 0, -1]) == lattice[2, 0, -2]
    @test lattice[1, 0, 0] == lattice[1, 0, 0]
    @test_broken lattice[1, 0, 0] == lattice[1, 0, 1]
    # @test_broken lattice[0, 1, 0] == make_grid_vector(expand(lattice), SVector(0, 1, 0))

    # Cartesian coordinates are already tested in `test_grid.jl`

    # The grid vectors of the dual grid.
    o_1 = [0, 1, 1] / sqrt(2)
    o_2 = [1, 0, 1] / sqrt(2)
    o_3 = [0, 0, 1] / sqrt(2)
    oblique_basis = (Vector3(o_1...), Vector3(o_2...), Vector3(o_3...))
    reciprocal_lattice = make_grid(ReciprocalLattice3D, oblique_basis, size_to_domain(sizes))
    homecell = transform_grid(reciprocal_lattice)
    @test isapprox(braket(dagger(homecell[0, 0, 1]), reciprocal_lattice[0, 0, 1]), π / 2, atol = 1e-7)

    ## Using a more convenient syntax.
    @test isapprox(homecell[0, 0, 1]' * reciprocal_lattice[0, 0, 1], π / 2, atol = 1e-7)

    ## Test that indexing a grid with the linear index of a grid vector gives back the grid vector.
    for r in homecell
        @test homecell(linear_index(r)) == r
    end
end

