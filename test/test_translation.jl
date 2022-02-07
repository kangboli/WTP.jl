using WTP
using Test


@testset "OnGrid" begin
    wave_functions_list = wave_functions_from_directory(joinpath(test_1_dir, "si.save"))
    ũ = orbital_set_from_save(wave_functions_list);
    brillouin_zone = grid(ũ)
    gamma_point = brillouin_zone[0, 0, 0]
    ψ = ũ[gamma_point][1]
    g = grid(ψ)

    # Indexing with a grid vector.
    # The indexing of `OnGrid` is not separately tested much because it's logic
    # is entirly in `miller_to_standard`, which is independently tested.
    @test ψ[g[0, 0, 0]] == elements(ψ)[1]
    
    # Translation.
    @test translation(ψ) == (0, 0, 0)
    translated = translate(ψ, [1, 1, 1])
    translation(translated) == (1, 1, 1)

    # Standard representation of translated orbitals.
    @time translated_std = standardize(translated);
    @test translation(translated_std) == (0, 0, 0)

    # Verify the translation.
    @test isapprox(translated_std[g[1, 1, 1]] - ψ[g[0, 0, 0]], 0, atol=1e-7)
    @test isapprox(translated_std[g[1, 0, 1]] - ψ[g[0, -1, 0]], 0, atol=1e-7)

    # Kpoint continuity test
    z2m1 = norm(braket(dagger(ũ[brillouin_zone[0, 0, 2]][1]),
    standardize(translate(ũ[brillouin_zone[0, 0, -1]][1], [0, 0, -1]))))
    z21 = norm(braket(dagger(ũ[brillouin_zone[0, 0, 2]][1]), ũ[brillouin_zone[0, 0, 1]][1]))
    @test isapprox(z2m1, z21, atol=1e-5)


end