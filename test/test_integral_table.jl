@testset "integral table" begin

    wave_functions_list = wave_functions_from_directory(joinpath(test_1_dir, "si.save"))
    wannier = wannier_from_save(wave_functions_list)
    g = grid(wannier)
    gamma_point = g[0, 0, 0]

    ## Find shells

    neighbor_shells = find_shells(wannier, 1)
    @test g[-1, -1, -1] in shells[1]
    @test g[1, 1, 1] in shells[1]
    @test g[0, 0, 1] in shells[1]
    @test g[0, -1, 0] in shells[1]

    ## Find weights

    @test isapprox(compute_weights(neighbor_shells)[1],  5.3360380375, atol=1e-6)

    scheme = W90FiniteDifference(wannier)
    integrals_from_wave_functions = neighbor_basis_integral(scheme)

    @test braket(dagger(wannier[gamma_point][1]), wannier[g[0, 0, 1]][1]) ==
          integrals_from_wave_functions[gamma_point, g[0, 0, 1]][1, 1]

    # integrals_from_wave_functions[g[0, 0, 1], gamma_point]
    # integrals(integrals_from_wave_functions)[g[0, 0, 1]=>gamma_point]
    # for (k1, k2) in keys(integrals(integrals_from_wave_functions))
    #     if k1 == gamma_point 
    #         println(k2)
    #     end
    # end


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


    # Check that this matches with the mmn file.

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

    amn = AMN(joinpath(test_1_dir, "output/pw2wan/si.amn"))
    gauge(amn, k_map)

end



# integral_from_mmn[g[-1,1,1], g[-2,1,1]][1, 2] - integrals_from_wave_functions[g[-1,1,1], g[-2,1,1]][1, 2]
# integral_from_mmn[g[-1,1,1], g[-2,1,1]][2, 2] - integrals_from_wave_functions[g[-1,1,1], g[-2,1,1]][6, 6]

# integral_from_mmn[g[-1,1,1], g[-2,1,1]][1, 1]
# integrals_from_wave_functions[g[-1,1,1], g[-2,1,1]][5, 5]

# diff =integral_from_mmn[g[-1,1,1], g[-2,1,1]] - integrals_from_wave_functions[g[-1,1,1], g[-2,1,1]][5:end, 5:end]
# diff[:,1]
# norm(integral_from_mmn[g[-1,1,1], g[-2,1,1]][1, 1] - integrals_from_wave_functions[g[-1,1,1], g[-2,1,1]][5, 5])

# wannier[g[-1,1,1]][5]
# braket(dagger(wannier[5, g[-1,1,1]]), wannier[5, g[-2,1,1]])

# braket(dagger(wannier[g[-1,1,1]][5]), 
# standardize(translate(wannier[g[-2,1,1]][5], reciprocal_lattice(wannier)[1, 0, 0])))
# j

# braket(dagger(wannier[g[-1,1,1]][5]), 
# standardize(translate(wannier[reset_overflow(g[-2,1,1])][5], 
# reciprocal_lattice(wannier)[1, 0, 0])))

# wannier[5, g[-2,1,1]]

# standardize(translate(wannier[g[-2,1,1]][5], reciprocal_lattice(wannier)[1, 0, 0]))
