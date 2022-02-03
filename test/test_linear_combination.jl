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
b_1 = DefaultLinearCombination([1, 0, 0], oblique_basis, true)
b_2 = DefaultLinearCombination([0, 1, 0], oblique_basis, true)
@test isapprox(basis_transform(coefficients(b_1), oblique_basis, CARTESIAN_BASIS), o_1, atol=1e-7 )
@test isapprox(basis_transform(coefficients(b_2), oblique_basis, CARTESIAN_BASIS), o_2, atol=1e-7 )


@testset "Linear Combinationination of Orbitals" begin

    wave_functions_list = wave_functions_from_directory(joinpath(test_1_dir, "si.save"))
    wannier = wannier_from_save(wave_functions_list);
    gamma_point = grid(wannier)[0, 0, 0]
    u1gamma = wannier[gamma_point][1]
    u2gamma = wannier[gamma_point][2]

    comb_1 = UnkOrbital(u1gamma)
    comb_2 = UnkOrbital(u2gamma)
    
    @test basis(comb_1)[1] == u1gamma
    @test coefficients(comb_1) == [1]

    comb_3 = add(comb_1, comb_2);

    @test coefficients(comb_3) == [1., 1.]
    @test coefficients(negate(comb_3))  == -[1., 1.]
    @test coefficients(add(comb_3, comb_3)) == [2., 2.]
    @test coefficients(minus(comb_3, comb_1)) == [0., 1.]
    @test coefficients(mul(2, comb_3)) == [2., 2.]

    @test isapprox(braket(dagger(comb_1), comb_1), 1, atol=1e-7)
    @test isapprox(braket(dagger(comb_2), comb_1), 0, atol=1e-7)

    ut1gamma = wannier[1, gamma_point]
    ut1kappa = wannier[1, grid(wannier)[0, 0, 1]]
    c1 = UnkOrbital(ut1gamma)
    c2 = UnkOrbital(ut1kappa)
    c3 = add(c1, c2, false)
    @test norm(braket(dagger(c3), c3)) < 4
end