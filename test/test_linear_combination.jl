"""
using Test
using WTP
"""

using LinearAlgebra

struct DefaultLinearCombination <: LinearCombination
    _coefficients
    basis
    ket::Bool
end

o_1 = [0, 1, 1] / sqrt(2)
o_2 = [1, 0, 1] / sqrt(2)
o3 = [0, 0, 1] / sqrt(2)
oblique_basis = (Vector3(o_1...), Vector3(o_2...), Vector3(o3...))
b1 = DefaultLinearCombination([1, 0, 0], oblique_basis, true)
b2 = DefaultLinearCombination([0, 1, 0], oblique_basis, true)
@test isapprox(coefficients(basis_transform(b1, CARTESIAN_BASIS)), o_1, atol=1e-7 )
@test isapprox(coefficients(basis_transform(b2, CARTESIAN_BASIS)), o_2, atol=1e-7 )


@testset "Linear Combinationination of Orbitals" begin
    wannier = wannier_from_save(joinpath(experiment_dir, "si.save"));
    gamma_point = grid(wannier)[0, 0, 0]
    u1gamma = wannier[gamma_point][1]
    u2gamma = wannier[gamma_point][2]

    comb1 = UnkOrbital(u1gamma)
    comb2 = UnkOrbital(u2gamma)
    
    @test basis(comb1)[1] == u1gamma
    @test coefficients(comb1) == [1]

    comb3 = add(comb1, comb2);

    @test coefficients(comb3) == [1., 1.]
    @test coefficients(negate(comb3))  == -[1., 1.]
    @test coefficients(add(comb3, comb3)) == [2., 2.]
    @test coefficients(minus(comb3, comb1)) == [0., 1.]
    @test coefficients(mul(2, comb3)) == [2., 2.]

    @test isapprox(braket(dagger(comb1), comb1), 1, atol=1e-7)
    @test isapprox(braket(dagger(comb2), comb1), 0, atol=1e-7)

    ut1gamma = wannier[1, gamma_point]
    ut1kappa = wannier[1, grid(wannier)[0, 0, 1]]
    c1 = UnkOrbital(ut1gamma)
    c2 = UnkOrbital(ut1kappa)
    c3 = add(c1, c2, false)
    @test norm(braket(dagger(c3), c3)) < 4
end