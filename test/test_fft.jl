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
        HomeCell(CARTESIAN_BASIS, size_to_domain((unk.nx, unk.ny, unk.nz))), gamma_point, 1)
    # wtp_normalize!(real_orbital)
    real_transformed = fft(real_orbital)
    norm(elements(real_orbital))
    norm(elements(real_transformed))
    # wtp_normalize!(real_transformed)
    transformed_back = ifft(real_transformed)
    norm(elements(transformed_back))
    @test norm(elements(transformed_back) - elements(real_orbital)) < 1e-7
    @test norm(elements(real_transformed) - elements(wannier[gamma_point][1]) ) < 1e-7
    

    # for n = 1:20
    #     real_transformed = ifft(real_orbitals[n])
    #     wtp_normalize!(real_transformed)
    #     @test (norm(elements(real_transformed) - elements(wannier[gamma_point][n]))) < 1e-7
    # end
end
