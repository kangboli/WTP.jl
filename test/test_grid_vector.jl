"""
using Test
using WTP
"""

reciprocal_basis = CARTESIAN_BASIS

sizes = (4, 4, 4)
lattice = ReciprocalLattice(reciprocal_basis, size_to_domain(sizes))
g = GridVector{ReciprocalLattice}(lattice, [0, 1, 2], true)
# A grid vector is defined on a grid.
@test grid(g) == lattice
# THe size of the grid vector should be that of the grid.
@test size(grid(g)) == sizes
# Indexing a grid should yield a grid vector.
@test lattice[0, 1, 2] == g
# Basic vector algebra.
@test add(lattice[0, 0, 1], lattice[0, 1, 0]) == lattice[0, 1, 1]

# The _coefficients should wrap around once they are out of the grid.
h = lattice[0, 0, 3]
@test wrapped(h) == [0, 0, -1]
# The third coefficient has been wrapped once.
@test overflow(h) == [0, 0, 1]
# Wrapping should compute with negation.
@test has_overflow(h)

# Get the image within the grid domain.
h2 = reset_overflow(h)
@test !has_overflow(h2)

# Get back the image out of the grid domain.
@test h == add_overflow(h2, overflow(h))

@test 9 == braket(dagger(h), h)

@test coefficients(negate(h)) == -(coefficients(h))

# Snapping after taking the coordinates should recover the vector
@test snap(lattice, [0, 1, 2]) == g

@testset "Miller indices" begin
    sizes = (4, 4, 4)
    @test [1, 1, 1] == miller_to_standard(sizes, [0, 0, 0], [0, 0, 0])
    @test [3, 4, 2] == miller_to_standard(sizes, [1, -1, 1], [1, 0, 0])

    @test [0, 0, 0] == standard_to_miller(sizes, [1, 1, 1], [0, 0, 0])
    @test [1, -1, 1] == standard_to_miller(sizes, [3, 4, 2], [1, 0, 0])

    one_dim_index = three_to_one(miller_to_standard(sizes, [1, -1, 1], [1, 0, 0])..., sizes)
    @test [1, -1, 1] == standard_to_miller(sizes, one_to_three(one_dim_index, sizes), [1, 0, 0])
end

# The grid vectors of the dual grid.
o_1 = [0, 1, 1] / sqrt(2)
o_2 = [1, 0, 1] / sqrt(2)
o_3 = [0, 0, 1] / sqrt(2)
oblique_basis = (Vector3(o_1...), Vector3(o_2...), Vector3(o_3...))
reciprocal_lattice = ReciprocalLattice(oblique_basis, size_to_domain(sizes))
homecell = transform_grid(reciprocal_lattice)
@test isapprox(braket(dagger(homecell[0,0,1]), reciprocal_lattice[0,0,1]),  π/2, atol=1e-7)

## Using a more convenient syntax.
@test isapprox(homecell[0,0,1]' * reciprocal_lattice[0,0,1],  π/2, atol=1e-7)

## Test that indexing a grid with the linear index of a grid vector gives back the grid vector.

for r in homecell
    @test homecell[linear_index(r)] == r
end


