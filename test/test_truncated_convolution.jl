## Load the wave function.

wave_functions_list = wave_functions_from_directory(joinpath(test_6_dir, "si.save"))
ũ = wannier_from_save(wave_functions_list);

brillouin_zone = grid(ũ)

## Get the map of kpoints.

k_map, _ = i_kpoint_map(wave_functions_list)


## Load the gauge transform amn and verify it.

amn = AMN(joinpath(test_6_dir, "output/pw2wan/si.amn"))
U = Gauge(grid(ũ), amn, k_map)

## Construct the scheme.

scheme = W90FiniteDifference3D(ũ)
G_1 = gauge_gradient(U, scheme, brillouin_zone, TruncatedConvolution) 
G_2 = gauge_gradient(U, scheme, brillouin_zone, BranchNaive)

sum(brillouin_zone) do k
    norm(G_1[k] - G_2[k])
end

sum(brillouin_zone) do k
    norm(G_1[k])
end

sum(brillouin_zone) do k
    norm(G_2[k])
end

## Apply the gauge transform on the basis integrals to get the neighbor integrals.

M = gauge_transform(neighbor_basis_integral(scheme), U)

## Get the center from the neighbor integrals and the scheme. Test it against the output 
# from Wannier90.

center(M, scheme, 1, TruncatedConvolution)
center(M, scheme, 2, TruncatedConvolution)
center(M, scheme, 3, TruncatedConvolution)
center(M, scheme, 4, TruncatedConvolution)

spread(M, scheme, 1, TruncatedConvolution)
spread(M, scheme, 1, BranchNaive)
spread(M, scheme, 2, TruncatedConvolution)
spread(M, scheme, 3, TruncatedConvolution)
spread(M, scheme, 4, TruncatedConvolution)

u = ifft(ũ)
u = set_gauge(u, U)
phase = phase_factors(u)

wannier_orbitals = u(:, phase)

supercell = expand(orbital_grid(u), [size(brillouin_zone)...])
r2 = map(supercell) do r
    norm(r)^2
end
r̃2 = fft(r2, false)
densities = abs2.(wannier_orbitals)

center_spread(fft(densities[1], false), r̃2)
center_spread(fft(densities[2], false), r̃2)
center_spread(fft(densities[3], false), r̃2)
center_spread(fft(densities[4], false), r̃2)