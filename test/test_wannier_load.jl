@testset "WannierSet" begin
    wave_functions_list = wave_functions_from_directory(joinpath(test_1_dir, "si.save"))
    wannier = wannier_from_save(wave_functions_list)
    g = grid(wannier)

    gamma_point = g[0,0,0]
    u1gamma = wannier[gamma_point][1]
    @test kpoint(u1gamma) == gamma_point
    @test index_band(u1gamma) == 1
    # Orthogonality test.
    @test isapprox(braket(dagger(wannier[gamma_point][1]), wannier[gamma_point][2]), 0, atol=1e-7)
    @test isapprox(braket(dagger(wannier[gamma_point][8]), wannier[gamma_point][16]), 0, atol=1e-7)

    # Uniqueness test.
    for k in g
        if isapprox(norm(braket(dagger(wannier[gamma_point][1]), wannier[k][1])), 1.0, atol=1e-3)
            @test k == gamma_point
        end
    end

    # Affinity test. Uₙₖ should be similar for Nearby k-points.
    z01 = norm(braket(dagger(wannier[gamma_point][1]), wannier[g[0, 0, 1]][1])) 
    zm10 = norm(braket(dagger(wannier[g[0, 0, -1]][1]), wannier[gamma_point][1]))
    @test z01 > 0.9

    z12 = norm(braket(dagger(wannier[g[0, 0, 1]][1]), wannier[g[0, 0, 2]][1]))
    
    @test z12 > 0.8

    # Gamma symmetry test. Uₙₖ should be symmetric around the gamma point.
    @test isapprox(z01, zm10, atol=1e-5)

    # Test reading from UNK files. Applying FFT to this should yield the same
    # result as reading from the wave_functions files.
    wannier2 = wannier_from_unk_dir(joinpath(test_1_dir, "unk"), wave_functions_list)

    for kpoint in g
        for b in 1:20
            @test isapprox(norm(elements(wannier2[kpoint][1]) - elements(wannier[kpoint][1])), 0, atol=1e-6)
        end
    end
    
end


    # wannier_e01 = wannier_from_save(joinpath(e01_dir, "si.save"))
    # wannier_e01_unk = wannier_from_unk_dir(joinpath(e01_dir, "unk"),
    #                                         joinpath(e01_dir, "si.save"))

    # g_e01 = grid(wannier_e01)
    # gamma_point = g_e01[0,0,0]
    # norm(elements(wannier_e01[gamma_point][1]) - elements(wannier_e01_unk[gamma_point][5]))