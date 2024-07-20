using WTP
using Test
using LinearAlgebra


@testset "fft" begin
    wave_functions_list = wave_functions_from_directory(joinpath(test_1_dir, "si.save"))
    ũ = orbital_set_from_save(wave_functions_list);
    gamma_point = grid(ũ)[0, 0, 0]

    # Load the real orbital from a unk file.
    unk = UNK(joinpath(test_1_dir, "output/pw2wan/unk/UNK00001.1"))
    real_orbital = single_orbital_from_unk(unk,
        make_grid(HomeCell3D, CARTESIAN_BASIS, size_to_domain((unk.nx, unk.ny, unk.nz))), gamma_point, 1)

    # Transform the real orbital to the reciprocal space and match 
    # with one read from the .save folder.
    real_transformed = fft(real_orbital)
    @test isapprox(norm(elements(real_orbital)), 1, atol=1e-12)
    @test isapprox(norm(elements(real_transformed)), 1, atol=1e-12)
    @test norm(elements(real_transformed) - elements(ũ[gamma_point][1]) ) < 1e-6

    # Transform it back and make sure that it matches the original.
    transformed_back = ifft(real_transformed)
    @test norm(elements(transformed_back) - elements(real_orbital)) < 1e-6

    # Test that the in-place fft mutates the orbital.
    real_transformed = @fft!(real_orbital)
    @test real_orbital === nothing
    # Test that the in-place ifft mutate the orbital back.
    real_orbital = @ifft!(real_transformed)
    @test real_transformed === nothing

    @ifft!(ũ)
    @test ũ === nothing
end
