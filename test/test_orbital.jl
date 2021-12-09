"""
using Test
using WTP
"""

@testset "Translation" begin
    wave_functions_list = wave_functions_from_directory(joinpath(test_1_dir, "si.save"))
    wannier = wannier_from_save(wave_functions_list);
    brillouin_zone = grid(wannier)
    gamma_point = brillouin_zone[0, 0, 0]
    u1gamma = wannier[gamma_point][1]
    g = grid(u1gamma)
    @test translation(u1gamma) == g[0, 0, 0]
    # Lazy translation.
    translated = translate(u1gamma, g[1, 1, 1])
    @test translation(translated) == g[1, 1, 1]

    # Convert the translation to standard representation.
    @time translated_std = standardize(translated);
    @test translation(translated_std) == g[0, 0, 0]

    # Verify the translation.
    @test isapprox(translated_std[g[1, 1, 1]] - u1gamma[g[0, 0, 0]], 0, atol=1e-7)
    @test isapprox(translated_std[g[1, 0, 1]] - u1gamma[g[0, -1, 0]], 0, atol=1e-7)
    # Kpoint continuity test
    z2m1 = norm(braket(dagger(wannier[brillouin_zone[0, 0, 2]][1]),
    standardize(translate(wannier[brillouin_zone[0, 0, -1]][1], g[0, 0, -1]))))

    z21 = norm(braket(dagger(wannier[brillouin_zone[0, 0, 2]][1]), wannier[brillouin_zone[0, 0, 1]][1]))
    @test isapprox(z2m1, z21, atol=1e-5)

end