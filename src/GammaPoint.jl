
# export phase_expectation

# """

# """
# abstract type GammaEvaluator end

# struct W90GammaEvaluator <: GammaEvaluator
#     x::Matrix{ComplexFxx}
#     y::Matrix{ComplexFxx}
#     z::Matrix{ComplexFxx}
#     cell_length::Number
# end

# """
# x = ⟨wₘ|exp(-i(2π/L)x)|wₙ⟩
# """
# function W90GammaEvaluator(u::Wannier{UnkBasisOrbital{ReciprocalLattice}})
#     gamma = grid(u)[0, 0, 0]
#     u = u[gamma]
#     homecell = transform_grid(grid(u[1]))
#     cell_length = (vector3_to_matrix(basis(homecell))*diagm([size(homecell)...]))[1, 1]

#     n_band = length(u)
#     x, y, z = (zeros(ComplexF64, n_band, n_band) for _ = 1:3)
#     for i = 1:n_band, j = 1:n_band
#         x[i, j] = u[i]' * (u[j] >> [-1, 0, 0])
#         y[i, j] = u[i]' * (u[j] >> [0, -1, 0])
#         z[i, j] = u[i]' * (u[j] >> [0, 0, -1])
#     end

#     W90GammaEvaluator(x, y, z, cell_length)
# end

# x(e::W90GammaEvaluator) = e.x
# y(e::W90GammaEvaluator) = e.y
# z(e::W90GammaEvaluator) = e.z
# cell_length(e::W90GammaEvaluator) = e.cell_length

# function gauge_transform(
#     x::Matrix{ComplexFxx},
#     y::Matrix{ComplexFxx},
#     z::Matrix{ComplexFxx},
#     u::Matrix{ComplexFxx},
# )
#     (adjoint(u) * x * u, adjoint(u) * y * u, adjoint(u) * z * u)
# end

# function center(
#     x::Matrix{ComplexFxx},
#     y::Matrix{ComplexFxx},
#     z::Matrix{ComplexFxx},
#     evaluator::W90GammaEvaluator,
#     n::Int,
# )
#     l = cell_length(evaluator)
#     c_x = -l / (2 * π) * imag(log(x[n, n]))
#     c_y = -l / (2 * π) * imag(log(y[n, n]))
#     c_z = -l / (2 * π) * imag(log(z[n, n]))
#     return (c_x, c_y, c_z)
# end

# function second_moment(
#     x::Matrix{ComplexFxx},
#     y::Matrix{ComplexFxx},
#     z::Matrix{ComplexFxx},
#     evaluator::W90GammaEvaluator,
#     n::Int,
# )
#     sum((m -> ((1 - abs2(m[n, n])) + imag(log(m[n, n]))^2)).([x, y, z]))
# end

# # """
# # ⟨wₘ|exp(-i(2π/L)x)|wₙ⟩
# # """
# # function phase_expectation(o_1::UnkBasisOrbital{ReciprocalLattice}, o_2::UnkBasisOrbital{ReciprocalLattice})
# #     reciprocal_lattice = grid(o_1)
# #     o_x = standardize(translate(o_2, reciprocal_lattice[-1, 0, 0]))
# #     o_y = standardize(translate(o_2, reciprocal_lattice[0, -1, 0]))
# #     o_z = standardize(translate(o_2, reciprocal_lattice[0, 0, -1]))
# #     return [braket(dagger(o_1), o) for o in [o_x, o_y, o_z]]
# # end