
## Load the wave function.

wave_functions_list = wave_functions_from_directory(joinpath(test_2_dir, "si.save"))
wannier = wannier_from_save(wave_functions_list);

brillouin_zone = grid(wannier)

## Get the map of kpoints.

k_map, _ = i_kpoint_map(wave_functions_list)


## Load the gauge transform amn and verify it.

amn = AMN(joinpath(test_2_dir, "output/pw2wan/si.amn"))
U = Gauge(grid(wannier), amn, k_map)

## Construct the scheme.

scheme = W90FiniteDifference(wannier)

## Apply the gauge transform on the basis integrals to get the neighbor integrals.

M = gauge_transform(neighbor_basis_integral(scheme), U)

## Get the center from the neighbor integrals and the scheme. Test it against the output 
# from Wannier90.

@test isapprox(center(M, scheme, 1), [1.751302, -3.656156, -3.154520], atol = 1e-6)
@test isapprox(center(M, scheme, 2), [0.662933, 1.224859, 0.588340], atol = 1e-6)
@test isapprox(center(M, scheme, 3), [0.751350, -1.252297, 0.334295], atol = 1e-6)
@test isapprox(center(M, scheme, 4), [0.745342, 0.390014, -1.239850], atol = 1e-6)
@test isapprox(center(M, scheme, 5), [-1.204118, 0.700036, 0.596825], atol = 1e-6)
@test isapprox(center(M, scheme, 6), [-0.731020, -0.808056, -1.032846], atol = 1e-6)
@test isapprox(center(M, scheme, 7), [2.542534, -2.518040, -1.119903], atol = 1e-6)
@test isapprox(center(M, scheme, 8), [2.213292, -1.116839, -2.496139], atol = 1e-6)
@test isapprox(center(M, scheme, 9), [-1.433266, -1.527580, 1.371041], atol = 1e-6)
@test isapprox(center(M, scheme, 10), [2.083056, -2.221819, -3.950472], atol = 1e-6)
@test isapprox(center(M, scheme, 11), [3.593463, -3.384827, -3.215746], atol = 1e-6)
@test isapprox(center(M, scheme, 12), [0.129847, -0.030452, 1.388568], atol = 1e-6)
@test isapprox(center(M, scheme, 13), [-1.283663, 1.192107, -1.304339], atol = 1e-6)
@test isapprox(center(M, scheme, 14), [1.120100, -2.677675, -2.025788], atol = 1e-6)
@test isapprox(center(M, scheme, 15), [3.978456, -2.243099, -2.060208], atol = 1e-6)
@test isapprox(center(M, scheme, 16), [-3.917846, 1.254606, 1.095461], atol = 1e-6)
@test isapprox(center(M, scheme, 17), [-1.054748, 4.053913, 1.282466], atol = 1e-6)
@test isapprox(center(M, scheme, 18), [-1.549755, 1.039485, 4.050678], atol = 1e-6)
@test isapprox(center(M, scheme, 19), [-3.226431, 4.489533, 3.791046], atol = 1e-6)
@test isapprox(center(M, scheme, 20), [-2.434481, 2.608875, 2.188386], atol = 1e-6)

## Get the spread from the neighbor integrals and the scheme. Test it against the output from
# Wannier90.

spread(n) = second_moment(M, scheme, n) - norm(center(M, scheme, n))^2

@test isapprox(spread(1), 15.47179726, atol = 1e-5)
@test isapprox(spread(2), 13.17038702, atol = 1e-5)
@test isapprox(spread(3), 7.64407335, atol = 1e-5)
@test isapprox(spread(4), 8.25678268, atol = 1e-5)
@test isapprox(spread(5), 9.77938941, atol = 1e-5)
@test isapprox(spread(6), 7.87085298, atol = 1e-5)
@test isapprox(spread(7), 9.17499993, atol = 1e-5)
@test isapprox(spread(8), 8.44416108, atol = 1e-5)
@test isapprox(spread(9), 9.83532580, atol = 1e-5)
@test isapprox(spread(10), 10.35456759, atol = 1e-5)
@test isapprox(spread(11), 11.23672086, atol = 1e-5)
@test isapprox(spread(12), 9.93107295, atol = 1e-5)
@test isapprox(spread(13), 8.86266731, atol = 1e-5)
@test isapprox(spread(14), 7.67576041, atol = 1e-5)
@test isapprox(spread(15), 10.03275833, atol = 1e-5)
@test isapprox(spread(16), 15.04579948, atol = 1e-5)
@test isapprox(spread(17), 14.16030120, atol = 1e-5)
@test isapprox(spread(18), 12.06222393, atol = 1e-5)
@test isapprox(spread(19), 13.83179659, atol = 1e-5)
@test isapprox(spread(20), 23.16652774, atol = 1e-5)























