
using WTP
using Test
using LinearAlgebra

@testset "W90 Gradient Descent" begin
    wave_functions_list = wave_functions_from_directory(joinpath(test_5_dir, "si.save"))
    ũ = orbital_set_from_save(wave_functions_list)
    brillouin_zone = grid(ũ)

    ## Get the map of kpoints.

    k_map, _ = i_kpoint_map(wave_functions_list)

    ## Load the gauge transform amn and verify it.
    amn = AMN(joinpath(test_5_dir, "output/pw2wan/si.amn"))
    U = Gauge(grid(ũ), amn, k_map)

    scheme = CosScheme3D(ũ)

    optimizer = ILAOptimizer(scheme)
    supercell = expand(transform_grid(orbital_grid(ũ)), [size(brillouin_zone)...])
    r̃2 = fft(map(r->norm(r)^2, supercell) , false)

    optimizer.meta[:r̃2] = r̃2
    optimizer.meta[:ũ] = ũ
    optimizer.meta[:ila_spread] = Vector{Vector{Float64}}()
    optimizer.meta[:convolutional_spread] = Vector{Vector{Float64}}()
    optimizer.meta[:center_difference] = Vector{Vector{Float64}}()
    U_optimal = optimizer(U, W90BranchCut; α_0=1)

    M_optimal = gauge_transform(neighbor_basis_integral(scheme),  U_optimal)
    @test all(i->isapprox(spread(M_optimal, scheme, i), 5.712, atol=1e-3), 1:4)
end