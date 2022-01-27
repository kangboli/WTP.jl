
## Load the wave function.

wave_functions_list = wave_functions_from_directory(joinpath(test_5_dir, "si.save"))
wannier = wannier_from_save(wave_functions_list);

brillouin_zone = grid(wannier)

## Get the map of kpoints.

k_map, _ = i_kpoint_map(wave_functions_list)


## Load the gauge transform amn and verify it.

amn = AMN(joinpath(test_5_dir, "output/pw2wan/si.amn"))
U = Gauge(grid(wannier), amn, k_map)

## Construct the scheme.

scheme = W90FiniteDifference3D(wannier)

## Apply the gauge transform on the basis integrals to get the neighbor integrals.

M = gauge_transform(neighbor_basis_integral(scheme), U)
# M = neighbor_basis_integral(scheme)

## Get the center from the neighbor integrals and the scheme. Test it against the output 
# from Wannier90.


G = ila_gradient(M, scheme, brillouin_zone)
Γ = brillouin_zone[0, 0, 0]


isapprox(G[Γ][1, 1], 0.0000000000000000 + -2.89588732047005790E-002im, atol=1e-7)
isapprox(G[Γ][1, 2], -6.50366255112471980E-002 + 2.33213503851723916E-002im, atol=1e-7)
isapprox(G[Γ][1, 3], -2.31003951750815596E-002 + -7.90234635323724260E-003im, atol=1e-7)
isapprox(G[Γ][1, 4], -5.84584658047603467E-004 + 1.91681715993616622E-002im, atol=1e-7)
isapprox(G[Γ][1, 5], 0.16094737978528517 + -1.41772686830576294E-002im, atol=1e-7)
isapprox(G[Γ][1, 6], -2.63131640845353371E-002 + 2.13773913586990799E-002im, atol=1e-7)
isapprox(G[Γ][1, 7], -3.77002949064640241E-002 + -5.80868261451320533E-002im, atol=1e-7)
isapprox(G[Γ][1, 8], 2.33641180322331533E-002 + 1.01266951368890051E-002im, atol=1e-7)
isapprox(G[Γ][1, 9], -3.31757136636838465E-002 + 6.51379468564497882E-002im, atol=1e-7)
isapprox(G[Γ][1, 10], 2.64706582767112418E-002 + 3.87672090718196027E-002im, atol=1e-7)
isapprox(G[Γ][1, 11], -3.76325063985818928E-002 + 2.92837234732843061E-003im, atol=1e-7)
isapprox(G[Γ][1, 12], 1.75055989451118220E-003 + -4.29272421004532126E-002im, atol=1e-7)
isapprox(G[Γ][1, 13], -4.55602119221194257E-002 + 6.56171138723575000E-003im, atol=1e-7)
isapprox(G[Γ][1, 14], 7.98572473712325925E-002 + -2.11161958604519050E-002im, atol=1e-7)
isapprox(G[Γ][1, 15], 0.12507339764252975 + 4.04150779575447507E-003im, atol=1e-7)
isapprox(G[Γ][1, 16], -3.79210897186743525E-003 + 5.74604816219117468E-003im, atol=1e-7)
isapprox(G[Γ][1, 17], -6.50094076416881869E-003 + -2.42858667979146577E-002im, atol=1e-7)
isapprox(G[Γ][1, 18], -9.01619251990430984E-003 + -7.49463872737146826E-003im, atol=1e-7)
isapprox(G[Γ][1, 19], 1.57189623177843470E-002 + 2.17366479667984377E-002im, atol=1e-7)
isapprox(G[Γ][1, 20], -2.74432369706255758E-002 + 1.36648862450632139E-002im, atol=1e-7)

optimizer = ILAOptimizer(scheme)

u = ifft(wannier)
phase = phase_factors(u)
supercell = expand(orbital_grid(u), [size(brillouin_zone)...])
r2 = map(supercell) do r
    norm(r)^2
end
r̃2 = fft(r2, false)

optimizer.meta[:r̃2] = r̃2
optimizer.meta[:u] = u
optimizer.meta[:phase] = phase
optimizer.meta[:ila_spread] = Vector{Vector{Float64}}()
optimizer.meta[:convolutional_spread] = Vector{Vector{Float64}}()
optimizer.meta[:center_difference] = Vector{Vector{Float64}}()
optimizer(U);

hcat(optimizer.meta[:ila_spread]...)

u = set_gauge(u, U)
wannier_orbitals = u(:, phase);
densities = abs2.(wannier_orbitals)
σ_total = 0

# center_spread(fft(densities[1], false), r̃2)
for ρ in densities
    c, σ = center_spread(fft(ρ, false), r̃2)
    println(c)
    println(σ)
    σ_total += σ
end
println(σ_total)

spread(M, scheme, 1)
spread(M, scheme, 2)
center(M, scheme, 2)
spread(M, scheme, 3)
spread(M, scheme, 4)