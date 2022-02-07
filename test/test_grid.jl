using WTP
using Test
using LinearAlgebra
using Profile

@testset "Grid" begin
    ## Create non-orthogonal homecells with even/odd number of grid points in each direction.
    homecell_basis = (Vector3(0, sqrt(3)/2, 0), Vector3(1, 1/2, sqrt(3)/2), Vector3(0, 0, 1/2))
    sizes_even = (4, 6, 8)
    homecell_even = make_grid(HomeCell3D, homecell_basis, size_to_domain(sizes_even))
    sizes_odd = (1, 3, 5)
    homecell_odd = make_grid(HomeCell3D, homecell_basis, size_to_domain(sizes_odd))
    @test domain(homecell_odd) == ((-0, 0), (-1, 1), (-2, 2))

    # Number of dimensions
    @test n_dims(typeof(homecell_even)) == 3

    # Get and set the domain.
    @test domain(homecell_even) == ((-2, 1), (-3, 2), (-4, 3))
    @test homecell_even == set_domain(homecell_even, size_to_domain(sizes_even))

    # Grid expansion
    expanded_even = expand(homecell_even, [2, 3, 4])
    @test domain(expanded_even) == ((-4, 3), (-9, 8), (-16, 15))
    expanded_odd = expand(homecell_odd, [1, 3, 2])
    @test domain(expanded_odd) == ((-0, 0), (-4, 4), (-5, 4))

    # Min/Max grid points
    @test mins(expanded_even) == (-4, -9, -16)
    @test maxes(expanded_even) == (3, 8, 15)

    # The center should be the origin.
    @test center(homecell_even) == (0, 0, 0)
    @test center(homecell_odd) == (0, 0, 0)

    # Transform the grid.
    lattice = transform_grid(homecell_even)

    # The product of homecell and reciprocal basis vectors 
    # should multiply to 2 π. That is A * B = 2 π.
    lattice_basis_mat = vector3_to_matrix(basis(lattice))
    homecell_basis_mat = vector3_to_matrix(basis(homecell_even))
    @test isapprox(lattice_basis_mat' * homecell_basis_mat * diagm([sizes_even...]), I * 2 * pi)
    ## Transform a grid twice should yield the original grid.
    @test domain(homecell_even) == domain(transform_grid(lattice))
    isapprox(homecell_basis_mat, vector3_to_matrix(basis(transform_grid(lattice))))

    # Iterating over the grid should yield the right number of grid vectors.
    @test length([g for g in lattice]) == prod(sizes_even)
    # Grid indexing

    @test coefficients(homecell_odd[0, 0, 0]) == [0, 0, 0]
    @test coefficients.(homecell_odd[0, 0, -1:1]) == [[0, 0, -1], [0, 0, 0], [0, 0, 1]]
    @test size(homecell_even[0, -1:1, -2:2]) == (3, 5)
    @test size(homecell_even[-1:1, -1:1, -2:2]) == (3, 3, 5)

    # Linear indexing.
    # Should start from the origin for a centered grid.
    @test coefficients(homecell_odd(1)) == [0, 0, 0]
    @test length(homecell_odd(1:4)) == 4
    @test length(homecell_odd(:)) == length(collect(homecell_odd))
    # The index order should match the iteration order.
    for (i, r) in enumerate(homecell_odd)
        @test homecell_odd(i) == r
    end

    # Snapping a cartesian coordinate to a grid point.
    @test snap(homecell_even, [0, sqrt(3), 0]) |> coefficients == [2, 0, 0]
    @test snap(homecell_even, [0, sqrt(3), 1.49]) |> coefficients == [2, 0, 3]

    # Grids that are not centered.
    homecell_off_center = make_grid(HomeCell3D, homecell_basis, ((-1, 2), (-3, 1), (-3, 7)))
    @test center(homecell_off_center)  == (1, -1, 2)

    # transform should still work.
    @test size(homecell_off_center) == size(transform_grid(homecell_off_center))
    # The cartesian coordinates are moved.
    @test cartesian(homecell_off_center[0, 0, 0]) == cartesian(homecell_even[1, -1, 2]) 
end