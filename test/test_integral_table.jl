using WTP
using Test
using LinearAlgebra

@testset "integral table" begin

    wave_functions_list = wave_functions_from_directory(joinpath(test_1_dir, "si.save"))
    ũ = orbital_set_from_save(wave_functions_list)
    g = grid(ũ)
    gamma_point = g[0, 0, 0]

    ## Find shells

    neighbor_shells = find_shells(grid(ũ), 1)
    @test g[-1, -1, -1] in neighbor_shells[1]
    @test g[1, 1, 1] in neighbor_shells[1]
    @test g[0, 0, 1] in neighbor_shells[1]
    @test g[0, -1, 0] in neighbor_shells[1]

    ## Find weights

    @test isapprox(compute_weights(neighbor_shells)[1],  5.3360380375, atol=1e-6)

    scheme = CosScheme3D(ũ)
    integrals_from_wave_functions = neighbor_basis_integral(scheme)

    @test braket(dagger(ũ[gamma_point][1]), ũ[g[0, 0, 1]][1]) ==
          integrals_from_wave_functions[gamma_point, g[0, 0, 1]][1, 1]


    # Testing against the mmn files.
    # 1    2    0    0    0
    @test isapprox(
        integrals_from_wave_functions[gamma_point, g[0, 0, 1]][5, 5],
        0.311042037957 + 0.173491015254im,
        atol = 1e-7,
    )

    @test isapprox(
        integrals_from_wave_functions[gamma_point, g[0, 0, 1]][6, 5],
        -0.171934685919 - 0.424207587459im,
        atol = 1e-7,
    )

    # 1    4    0    0   -1
    @test isapprox(
        integrals_from_wave_functions[gamma_point, g[0, 0, -1]][5, 5],
        -0.297360912008 - 0.196017199352im,
        atol = 1e-7,
    )


    # Check that our integral matches with the mmn file.

    k_map, _ = i_kpoint_map(wave_functions_list)
    mmn = MMN(joinpath(test_1_dir, "output/pw2wan/si.mmn"))
    integral_from_mmn = NeighborIntegral(mmn, k_map)

    for (kpoint, neighbor) in keys(integrals(integral_from_mmn))
        @test isapprox(
            norm(
                integral_from_mmn[kpoint, neighbor] -
                integrals_from_wave_functions[kpoint, neighbor],
            ),
            0,
            atol = 1e-6,
        )
    end

    # amn = AMN(joinpath(test_1_dir, "output/pw2wan/si.amn"))
    # gauge(amn, k_map)

end
