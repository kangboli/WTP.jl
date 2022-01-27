## Load the wave function.

wave_functions_list = wave_functions_from_directory(joinpath(test_5_dir, "si.save"))
ũ = wannier_from_save(wave_functions_list);
u = ifft(ũ)
brillouin_zone = grid(ũ)

## Get the map of kpoints.

k_map, _ = i_kpoint_map(wave_functions_list)

phase = phase_factors(u)

supercell = expand(orbital_grid(u), [size(brillouin_zone)...])

## Load the gauge transform amn and verify it.
amn = AMN(joinpath(test_5_dir, "output/pw2wan/si.amn"))
U = Gauge(grid(u), amn, k_map)

u = set_gauge(u, U)
wannier_orbitals = u(:, phase)
ρ = abs2(wannier_orbitals[1])
ρ̃ = fft(ρ, false)
r_min, σ_min = center_spread(ρ̃, r̃2)

## Construct the scheme.

scheme = W90FiniteDifference3D(ũ)
## Apply the gauge transform on the basis integrals to get the neighbor integrals.

M = gauge_transform(neighbor_basis_integral(scheme), U)
