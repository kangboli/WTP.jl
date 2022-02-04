using WTP
using Printf
using Test
using LinearAlgebra

@testset "Benzene Gamma" begin

    wave_functions_list = wave_functions_from_directory(joinpath(test_3_dir, "benzene.save"))
    wannier = wannier_from_save(wave_functions_list)
    gamma_point = grid(wannier)[0, 0, 0]

    k_map, _ = i_kpoint_map(wave_functions_list)

    ## Get the map of kpoints.

    amn = AMN(joinpath(test_3_dir, "output/pw2wan/benzene.amn"))
    U = Gauge(grid(wannier), amn, k_map)

    scheme = W90FiniteDifference3D(wannier, 1)

    M = gauge_transform(neighbor_basis_integral(scheme), U)

    @test isapprox(center(M, scheme, 1, BranchNaive), [13.664088, 13.584512, -14.820745], atol = 1e-6)
    @test isapprox(center(M, scheme, 2, BranchNaive), [12.299523, 13.400594, -14.930998], atol = 1e-6)
    @test isapprox(center(M, scheme, 3, BranchNaive), [11.440610, 12.624274, 14.740192], atol = 1e-6)
    @test isapprox(center(M, scheme, 4, BranchNaive), [13.392026, -13.454062, -14.966689], atol = 1e-6)
    @test isapprox(center(M, scheme, 5, BranchNaive), [13.090875, 13.075826, 14.574919], atol = 1e-6)
    @test isapprox(center(M, scheme, 6, BranchNaive), [12.541850, -13.527836, 14.994046], atol = 1e-6)
    @test isapprox(center(M, scheme, 7, BranchNaive), [14.258706, -14.931835, 14.836772], atol = 1e-6)
    @test isapprox(center(M, scheme, 8, BranchNaive), [11.183312, -13.922845, -14.944576], atol = 1e-6)
    @test isapprox(center(M, scheme, 9, BranchNaive), [12.393579, -14.493156, -14.765461], atol = 1e-6)
    @test isapprox(center(M, scheme, 10, BranchNaive), [12.635387, 13.903749, 14.884168], atol = 1e-6)
    @test isapprox(center(M, scheme, 11, BranchNaive), [-13.521396, 14.987840, -14.760481], atol = 1e-6)
    @test isapprox(center(M, scheme, 12, BranchNaive), [12.423043, -14.489749, 14.921353], atol = 1e-6)
    @test isapprox(center(M, scheme, 13, BranchNaive), [13.438646, -13.632424, 14.972581], atol = 1e-6)
    @test isapprox(center(M, scheme, 14, BranchNaive), [11.351283, 14.274368, -14.718336], atol = 1e-6)
    @test isapprox(center(M, scheme, 15, BranchNaive), [14.413077, -12.416998, 14.972585], atol = 1e-6)

    spread(n) = second_moment(M, scheme, n) - norm(center(M, scheme, n, BranchNaive))^2

    @test isapprox(spread(1), 11.08439110, atol = 1e-6)
    @test isapprox(spread(2), 5.85533247, atol = 1e-6)
    @test isapprox(spread(3), 7.79132358, atol = 1e-6)
    @test isapprox(spread(4), 7.32185066, atol = 1e-6)
    @test isapprox(spread(5), 6.61482690, atol = 1e-6)
    @test isapprox(spread(6), 7.72097925, atol = 1e-6)
    @test isapprox(spread(7), 10.58775234, atol = 1e-6)
    @test isapprox(spread(8), 10.71423090, atol = 1e-6)
    @test isapprox(spread(9), 11.71253115, atol = 1e-6)
    @test isapprox(spread(10), 10.32006113, atol = 1e-6)
    @test isapprox(spread(11), 4.82512068, atol = 1e-6)
    @test isapprox(spread(12), 7.28097787, atol = 1e-6)
    @test isapprox(spread(13), 8.11119769, atol = 1e-6)
    @test isapprox(spread(14), 5.67303440, atol = 1e-6)
    @test isapprox(spread(15), 6.77265437, atol = 1e-6)

end

# u_1 = wannier[gamma_point][1]
# reciprocal_lattice = grid(u_1)

# # phase_expectation(wannier[gamma_point][1], wannier[gamma_point][1])
# @test isapprox(abs(m_1 - (-0.84829899126445829-0.37973908322929639im)), 0, atol=1e-8)
# @test isapprox(abs(m_2 - (-0.92944040765360258+5.74048277679958017E-006im)), 0, atol=1e-8)
# @test isapprox(abs(m_3 - (-0.98288158837556283+6.10059947737300429E-007im)), 0, atol=1e-8)

# homecell = transform_grid(reciprocal_lattice)
# L = (vector3_to_matrix(basis(homecell)) * diagm([size(homecell)...]))[1,1]
# n_band = length(wannier[gamma_point])

# X = zeros(ComplexF64, n_band, n_band)
# Y = zeros(ComplexF64, n_band, n_band)
# Z = zeros(ComplexF64, n_band, n_band)

# for i = 1:n_band, j = 1:n_band
#     X[i, j], Y[i, j], Z[i, j] =  phase_expectation(wannier[gamma_point][i], wannier[gamma_point][j])
# end

# x_2 = adjoint(u[gamma_point]) * x * u[gamma_point]
# y_2 = adjoint(u[gamma_point]) * y * u[gamma_point]
# z_2 = adjoint(u[gamma_point]) * z * u[gamma_point]

# X_2 = adjoint(U[gamma_point]) * X_pw2wannier * U[gamma_point]
# Y_2 = adjoint(U[gamma_point]) * Y_pw2wannier * U[gamma_point]
# Z_2 = adjoint(U[gamma_point]) * Z_pw2wannier * U[gamma_point]


# for n = 1:n_band
#     c_x = - L / (2 * π) * imag(log(X_2[n, n]))
#     c_y = - L / (2 * π) * imag(log(Y_2[n, n]))
#     c_z = - L / (2 * π) * imag(log(Z_2[n, n]))

#     println(@sprintf("c_x: %.4f, ", c_x), @sprintf("c_y: %.4f, ", c_y), @sprintf("c_z: %.4f", c_z))
# end


# mmn = MMN(joinpath(test_3_dir, "output/pw2wan/benzene.mmn"), true)

# X_pw2wannier = mmn.integrals[1=>1]
# Y_pw2wannier = mmn.integrals[1=>2]
# Z_pw2wannier = mmn.integrals[1=>3]

# X_pw2wannier[1,1] 
# X[1,1] 


# U[gamma_point][1, 2]


# neighbor_integral = NeighborIntegral(mmn, k_map)

# # Perform the gauge transform and verify that the result is correct.
# neighbor_integral = gauge_transform(neighbor_integral, U)

# for r in homecell
#     u2_real[r] = exp.(-1im * (2 * π / L) * cartesian(r)[1]) * u2_real[r]
# end
# X = mmn.integrals[1=>1]
# Y = mmn.integrals[1=>2]
# Z = mmn.integrals[1=>3]
