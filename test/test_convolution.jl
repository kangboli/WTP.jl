using WTP
using Test
using LinearAlgebra
## Load the wave function.

wave_functions_list = wave_functions_from_directory(joinpath(test_6_dir, "si.save"))
ũ = wannier_from_save(wave_functions_list);
brillouin_zone = grid(ũ)

## Get the map of kpoints.

k_map, _ = i_kpoint_map(wave_functions_list)


## Load the gauge transform amn and verify it.
amn = AMN(joinpath(test_6_dir, "output/pw2wan/si.amn"))
U = Gauge(grid(ũ), amn, k_map)
ũ = set_gauge(ũ, U)

supercell = expand(transform_grid(orbital_grid(ũ)), [size(brillouin_zone)...])
r̃2 = fft(map(r->norm(r)^2, supercell), false)
norm(supercell[1, 0, 0])
supercell[1, 0, 0]
wannier_orbitals = commit_gauge(ũ)(:)
densities = abs2.(ifft.(wannier_orbitals))
sum(supercell) do r 
    densities[1][r] * norm(cartesian(r) + [1.18, 1.23, 1.23])^2
end
cartesian(supercell[1, 0, 0])
supercell
reciprocal_densities = (ρ->fft(ρ, false)).(densities)

## Construct the scheme.
scheme = W90FiniteDifference3D(ũ)
## Apply the gauge transform on the basis integrals to get the neighbor integrals.

M = gauge_transform(neighbor_basis_integral(scheme), U)

center_spread(reciprocal_densities[1], r̃2)
center_spread(reciprocal_densities[2], r̃2)
center_spread(reciprocal_densities[3], r̃2)
center_spread(reciprocal_densities[4], r̃2)
center(M, scheme, 1, BranchStable)
center(M, scheme, 2, BranchStable)
center(M, scheme, 3, BranchStable)
center(M, scheme, 4, BranchStable)

WTP.spread(M, scheme, 1)
WTP.spread(M, scheme, 2)
WTP.spread(M, scheme, 3)
WTP.spread(M, scheme, 4)
