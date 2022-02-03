"""
using Test
using WTP
"""

@testset "fft" begin
    wave_functions_list = wave_functions_from_directory(joinpath(test_1_dir, "si.save"))
    wannier = wannier_from_save(wave_functions_list);
    gamma_point = grid(wannier)[0, 0, 0]

    # Transform the UNK back to the reciprocal space.
    unk = UNK(joinpath(test_1_dir, "unk/UNK00001.1"))
    real_orbital = single_orbital_from_unk(unk,
        make_grid(HomeCell3D, CARTESIAN_BASIS, size_to_domain((unk.nx, unk.ny, unk.nz))), gamma_point, 1)
    # wtp_normalize!(real_orbital)
    real_transformed = fft(real_orbital)
    norm(elements(real_orbital))
    norm(elements(real_transformed))
    # wtp_normalize!(real_transformed)
    transformed_back = ifft(real_transformed)
    norm(elements(transformed_back))
    @test norm(elements(transformed_back) - elements(real_orbital)) < 1e-6
    @test norm(elements(real_transformed) - elements(wannier[gamma_point][1]) ) < 1e-6

    # Test that the in-place fft mutates the orbital.
    real_transformed = @fft!(real_orbital)
    @test real_orbital === nothing
    # Test that the in-place ifft mutate the orbital back.
    real_orbital = @ifft!(real_transformed)
    @test real_transformed === nothing

    @ifft!(wannier)
    @test wannier === nothing
end
