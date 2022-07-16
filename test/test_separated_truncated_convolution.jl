
using WTP
using Test
using LinearAlgebra
using SCDM
# using GLMakie
using UnicodePlots

wave_functions_list = wave_functions_from_directory(joinpath(test_6_dir, "si.save"))
ũ = orbital_set_from_save(wave_functions_list)
brillouin_zone = grid(ũ)
u = ifft(ũ)
u = commit_gauge(set_gauge(u, U_optimal))

homecell = orbital_grid(u)

## Get the map of kpoints.

k_map, _ = i_kpoint_map(wave_functions_list)

## Load the gauge transform amn and verify it.
amn = AMN(joinpath(test_6_dir, "output/pw2wan/si.amn"))
U = Gauge(grid(ũ), amn, k_map)
U = Gauge(grid(ũ), 4)

## Construct the scheme and compare the gradients.
scheme = CosScheme3D(ũ)

optimizer = ILAOptimizer(scheme)
U_optimal = optimizer(gauge(ũ), TruncatedConvolution, FletcherReeves; ϵ=1e-9, logging=true)

shells(scheme)
M = neighbor_basis_integral(scheme)
U_scdm, col = scdm_condense_phase(u, collect(1:4))
# M = gauge_transform(neighbor_basis_integral(scheme), U_scdm)
M = gauge_transform(neighbor_basis_integral(scheme), U_optimal)
Ψ_dagger = hcat(vectorize.(u[brillouin_zone[0, 0, 0]])...)
Ψ_dagger_next = hcat(vectorize.(u[brillouin_zone[0, 0, 1]])...)

phase = map(r->-1im*exp(1im * (brillouin_zone[0, 0, 1]' * r)), homecell)
# phase = map(r->exp(-1im * (brillouin_zone[0, 0, 2]' * r)), homecell)
-angle.(phase[homecell[0, :, :]]) + angle.(u[brillouin_zone[0, 0, 0]][1][homecell[0, :, :]]) - angle.(u[brillouin_zone[0, 0, 1]][1][homecell[0, :, :]]) |> UnicodePlots.heatmap
angle.(u[brillouin_zone[0, 0, 0]][1][homecell[0, :, :]])  |> UnicodePlots.heatmap
angle.(u[brillouin_zone[0, 0, -1]][1][homecell[0, :, :]])  |> UnicodePlots.heatmap
angle.(u[brillouin_zone[0, 0, 1]][1][homecell[0, :, :]])  |> UnicodePlots.heatmap
angle.(u[brillouin_zone[0, 0, -2]][1][homecell[0, :, :]])  |> UnicodePlots.heatmap
angle.(u[brillouin_zone[-1, 0, 0]][1][homecell[0, :, :]])  |> UnicodePlots.heatmap
angle.(phase[homecell[0, :, :]]) + angle.(u[brillouin_zone[0, 0, 1]][1][homecell[0, :, :]])  |> UnicodePlots.heatmap
angle.(phase[homecell[0, :, :]])  |> UnicodePlots.heatmap
angle.(u[brillouin_zone[0, 0, -2]][1][homecell[0, :, :]]) - angle.(phase[homecell[0, :, :]])  |> UnicodePlots.heatmap

circshift(angle.(u[brillouin_zone[0, 0, -1]][1][homecell[0, 0, :]]), 11 )  |> lineplot
angle.(u[brillouin_zone[0, 0, 1]][1][homecell[0, 0, :]])  |> lineplot
angle.(u[brillouin_zone[0, 0, -1]][1][homecell[0, 0, :]])  |> lineplot
angle.(u[brillouin_zone[0, 0, -2]][1][homecell[0, 0, :]])  |> lineplot
angle.(phase[homecell[0, 0, :]])  |> scatterplot
angle.(phase[homecell(col[1])]) 
phase[homecell(col[1])]
(angle.(u[brillouin_zone[0, 0, -1]][1][homecell[0, 0, :]])  - angle.(u[brillouin_zone[0, 0, -2]][1][homecell[0, 0, :]]))[1]
(angle.(u[brillouin_zone[1, -2, -1]][1][homecell[0, 0, :]])  - angle.(u[brillouin_zone[1, -2, -2]][1][homecell[0, 0, :]])) |> lineplot
phase = map(r->exp(1im * (brillouin_zone[0, 0, 1]' * r)), homecell)
println(norm(u[brillouin_zone[0, 0, -1]][1] * phase - u[brillouin_zone[0, 0, -2]][1] ))
println(norm(u[brillouin_zone[0, 0, -1]][1] * exp(1im * (-1.76)) - u[brillouin_zone[0, 0, -2]][1]))
for i = -pi:0.02:0
println("i: $(i)")
println(norm(u[brillouin_zone[0, 0, -1]][1] * exp(1im * i) - u[brillouin_zone[0, 0, -2]][1] ))
end
angle(phase[homecell[9, -3, -3]])
angle.(M11)
angle.(M12)
angle(phase[snap(homecell, center(M, scheme, 1, TruncatedConvolution)) ])

norm(elements(abs2(u[brillouin_zone[0, 0, -1]][1]) - abs2(u[brillouin_zone[0, 0, -2]][1])))
norm(elements(abs2(u[brillouin_zone[0, 0, -1]][1])))
1.9 + 4.22

angle(u[brillouin_zone[0, 0, -1]][1]' * u[brillouin_zone[0, 0, -2]][1])

abs.(u[brillouin_zone[1, 1, -1]][1][homecell[0, 0, :]])  |> lineplot
abs.(u[brillouin_zone[1, 1, -2]][1][homecell[0, 0, :]])  |> lineplot

angle.(u[brillouin_zone[-2, 0, -1]][1][homecell[0, 0, :]])  |> lineplot
angle.(u[brillouin_zone[-2, 0, -2]][1][homecell[0, 0, :]])  |> lineplot
abs.(u[brillouin_zone[-2, 1, -2]][1][homecell[0, 0, :]])  |> lineplot

angle.(u[brillouin_zone[0, 0, 0]][1][homecell[0, 0, :]])  |> lineplot
abs.(u[brillouin_zone[0, 0, 0]][1][homecell[-3, 9, :]])  |> lineplot
angle.(u[brillouin_zone[0, 0, 1]][1][homecell[0, 0, :]])  |> lineplot
homecell(col[1])
real.(u[brillouin_zone[0, 0, 1]][1][homecell[-3, 9, :]])  |> lineplot
imag.(u[brillouin_zone[0, 0, 1]][1][homecell[-3, 9, :]])  |> lineplot
u[brillouin_zone[0, 0, -2]][1]' * u[brillouin_zone[0, 0, -1]][1]
u[brillouin_zone[0, 0, 0]][1]' * u[brillouin_zone[0, 0, 1]][1]
M[brillouin_zone[0, 0, 1], brillouin_zone[0, 0, 2]]
norm(imag.(vectorize(phase * u[brillouin_zone[0, 0, 2]][1])))

real.(vectorize(u[brillouin_zone[0, 0, 0]][1]))' * real.(vectorize(u[brillouin_zone[0, 0, 1]][1]))


u[brillouin_zone[0, 0, 0]][1]' * 
norm(u[brillouin_zone[0, 0, 1]][1] - (u[brillouin_zone[0, 0, 0]][1] * phase))

abs.(u[brillouin_zone[0, 0, 2]][1][homecell[0, 0, :]])  |> lineplot
u[brillouin_zone[0, 0, 2]][1][homecell[0, 0, :]]

(angle.(Ψ_dagger) - angle.(Ψ_dagger_next))[:, 1] |> scatterplot

Ψ_dagger[col[1], :]  - Ψ_dagger_next[col[1], :] * exp(-1im * (brillouin_zone[0, 0, 1]' * homecell(col[1])))

Ψ_dagger[col[1], :]' * Ψ_dagger_next[col[1], :] 
angle(U_scdm[brillouin_zone[0, 0, 0]][1, :]' * U_scdm[brillouin_zone[0, 0, 1]][1, :])
angle.(M11)


U_scdm[brillouin_zone[0, 0, 0]]
U_scdm[brillouin_zone[0, 1, 0]]
U_scdm[brillouin_zone[0, 0, 1]]

M[brillouin_zone[0, 0, 0], brillouin_zone[0, 0, 1]]' * M[brillouin_zone[0, 0, 0], brillouin_zone[0, 0, 1]]
M11 = M[brillouin_zone[0, 0, 0], brillouin_zone[0, 0, 1]]
M12 = M[brillouin_zone[0, 0, 0], brillouin_zone[0, 0, -1]]
M21 = M[brillouin_zone[0, 0, 0], brillouin_zone[0, 1, 0]]
M22 = M[brillouin_zone[0, 0, 0], brillouin_zone[0, -1, 0]]
M31 = M[brillouin_zone[0, 0, 0], brillouin_zone[1, 0, 0]]
M32 = M[brillouin_zone[0, 0, 0], brillouin_zone[-1, 0, 0]]

R11 = hcat((n->M11[:, n] * M11[n, n]'/abs(M11[n, n]')).(1:4)...)
R12 = hcat((n->M12[:, n] * M12[n, n]'/abs(M12[n, n]')).(1:4)...)
R21 = hcat((n->M21[:, n] * M21[n, n]'/abs(M21[n, n]')).(1:4)...)
R22 = hcat((n->M22[:, n] * M22[n, n]'/abs(M22[n, n]')).(1:4)...)
R31 = hcat((n->M31[:, n] * M31[n, n]'/abs(M31[n, n]')).(1:4)...)
R32 = hcat((n->M32[:, n] * M32[n, n]'/abs(M32[n, n]')).(1:4)...)

R11 + R12
R21 + R22 - (R21 + R22)'
R31 + R32
R11 + R21 + R31 + R41 - (R12 + R22 + R32 + R42)
R11 - R12
R21 - R22
R11 - R11'
gauge_gradient_k_point_contribution(U, scheme, brillouin_zone, brillouin_zone[-1, 0, 0], STDC)
gauge_gradient_k_point_contribution(U, scheme, brillouin_zone, brillouin_zone[1, 0, 0], STDC)
gauge_gradient_k_point_contribution(U, scheme, brillouin_zone, brillouin_zone[0, 1, 0], STDC)
gauge_gradient_k_point_contribution(U, scheme, brillouin_zone, brillouin_zone[0, 0, 0], STDC)

M11 = M[brillouin_zone[-1, 0, 0], brillouin_zone[-1, 0, 1]]
M12 = M[brillouin_zone[-1, 0, 0], brillouin_zone[-1, 0, -1]]
M21 = M[brillouin_zone[-1, 0, 0], brillouin_zone[-1, 1, 0]]
M22 = M[brillouin_zone[-1, 0, 0], brillouin_zone[-1, -1, 0]]
M31 = M[brillouin_zone[-1, 0, 0], brillouin_zone[0, 0, 0]]
M32 = M[brillouin_zone[-1, 0, 0], brillouin_zone[-2, 0, 0]]
M41 = M[brillouin_zone[-1, 0, 0], brillouin_zone[0, 1, 1]]
M42 = M[brillouin_zone[-1, 0, 0], brillouin_zone[-2, -1, -1]]

R11 = hcat((n->M11[:, n] * M11[n, n]'/abs(M11[n, n]')).(1:4)...)
R12 = hcat((n->M12[:, n] * M12[n, n]'/abs(M12[n, n]')).(1:4)...)
R21 = hcat((n->M21[:, n] * M21[n, n]'/abs(M21[n, n]')).(1:4)...)
R22 = hcat((n->M22[:, n] * M22[n, n]'/abs(M22[n, n]')).(1:4)...)
R31 = hcat((n->M31[:, n] * M31[n, n]'/abs(M31[n, n]')).(1:4)...)
R32 = hcat((n->M32[:, n] * M32[n, n]'/abs(M32[n, n]')).(1:4)...)
R41 = hcat((n->M41[:, n] * M41[n, n]'/abs(M41[n, n]')).(1:4)...)
R42 = hcat((n->M42[:, n] * M42[n, n]'/abs(M42[n, n]')).(1:4)...)


R - R'
svd(M[brillouin_zone[0, 0, 0], brillouin_zone[0, 0, 1]])
svd(M[brillouin_zone[0, 0, 0], brillouin_zone[0, 1, 0]])
svd(M[brillouin_zone[0, 0, 0], brillouin_zone[1, 0, 0]])

svd(M[brillouin_zone[0, 0, -1], brillouin_zone[0, 0, 0]])
# @test sum(k -> norm(G_1[k] - G_2[k]), brillouin_zone) < 1
U_optimal, _ = scdm_condense_phase(ifft(ũ), collect(1:4))
U_optimal[brillouin_zone[0, 0, 0]] ./ U_optimal[brillouin_zone[0, 0, 1]]
U_optimal[brillouin_zone[0, 0, 0]]
U_optimal[brillouin_zone[0, 0, 1]] 
U_optimal[brillouin_zone[0, 0, 2]]
## Apply the gauge transform on the basis integrals to get the neighbor integrals.

M = gauge_transform(neighbor_basis_integral(scheme), U);
M = gauge_transform(neighbor_basis_integral(scheme), U_optimal);
M[brillouin_zone[0,0,0], brillouin_zone[0, 0, 1]]
M[brillouin_zone[0,0,-1], brillouin_zone[0, 0, 0]]
M_left = left_gauge_transform(neighbor_basis_integral(scheme), U_optimal, collect(brillouin_zone))
M_right = left_gauge_transform(neighbor_basis_integral(scheme), U_optimal, collect(brillouin_zone))
spread(M_left, M_right, scheme, 1, STDC)
spread(M, scheme, 2, SeparatedTruncatedConvolution)
spread(M, scheme, 3, SeparatedTruncatedConvolution)
spread(M, scheme, 4, SeparatedTruncatedConvolution)
spread(M, scheme, 1, TruncatedConvolution)
center(M, scheme, 1, TruncatedConvolution)
spread(M, scheme, 1, W90BranchCut)

## Get the center from the neighbor integrals and the scheme. Test it against the output 
# from Wannier90.

for i = 1:4
    c_1 = center(M, scheme, i, TruncatedConvolution)
    s_1 = spread(M, scheme, i, TruncatedConvolution)
    c_2 = center(M, scheme, i, W90BranchCut)
    s_2 = spread(M, scheme, i, W90BranchCut)
    @test isapprox(norm(c_1 - c_2), 0, atol = 1e-7)
    @test isapprox(s_1 - s_2, 0, atol = 0.5)
end