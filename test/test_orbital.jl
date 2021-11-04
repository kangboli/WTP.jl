"""
using Test
using WTP
"""


@testset "Translation" begin
    wannier = wannier_from_save(joinpath(experiment_dir, "si.save"));
    gamma_point = grid(wannier)[0, 0, 0]
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
    z2m1 = norm(braket(dagger(wannier[g[0, 0, 2]][1]),
    standardize(translate(wannier[g[0, 0, -1]][1], [0, 0, -1]))))

    z21 = norm(braket(dagger(wannier[g[0, 0, 2]][1]), wannier[g[0, 0, 1]][1]))
    @test isapprox(z2m1, z21, atol=1e-5)
end