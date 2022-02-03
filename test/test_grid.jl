
@testset "Grid" begin
    ## Create a non-orthogonal homecell.
    homecell_basis = (Vector3(0, sqrt(3)/2, 0), Vector3(1, 1/2, sqrt(3)/2), Vector3(0, 0, 1/2))
    sizes = (4, 6, 8)
    homecell = make_grid(HomeCell3D, homecell_basis, size_to_domain(sizes))
    lattice = transform_grid(homecell)

    # The product of homecell and reciprocal basis vectors 
    # should multiply to 2 π. That is A * B = 2 π.
    lattice_basis_mat = vector3_to_matrix(basis(lattice))
    homecell_basis_mat = vector3_to_matrix(basis(homecell))
    @test isapprox(lattice_basis_mat' * homecell_basis_mat * diagm([sizes...]), I * 2 * pi)

    ## Transform a grid twice should yield the original grid.
    @test domain(homecell) == domain(transform_grid(lattice))
    isapprox(homecell_basis_mat, vector3_to_matrix(basis(transform_grid(lattice))))

    # Iterating over the grid should yield the right number of grid vectors.
    @test length([g for g in lattice]) == prod(sizes)

    # Iterating should work for domains not centered as well.
    homecell_2 = make_grid(HomeCell3D, homecell_basis, ((-1, 2), (-3, 0), (-1, 2)))
    @test length([g for g in homecell_2]) == 64

    norms_of_g = map(g->norm(cartesian(g)), lattice)
    @test norms_of_g[lattice[0, 0, 0]] == 0
    
end