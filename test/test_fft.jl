"""
using Test
using WTP
"""

@testset "fft" begin
    wannier = wannier_from_save(joinpath(test_1_dir, "si.save"));
    gamma_point = grid(wannier)[0, 0, 0]


    # Transform the UNK back to the reciprocal space.
    unk = UNK(joinpath(test_1_dir, "unk/UNK00001.1"))
    real_orbitals = orbitals_from_unk(unk,
        HomeCell(CARTESIAN_BASIS, size_to_domain((unk.nx, unk.ny, unk.nz))), gamma_point)

    
    for n = 1:20
        real_transformed = fft(real_orbitals[n])
        wtp_normalize!(real_transformed)
        @test (norm(elements(real_transformed) - elements(wannier[gamma_point][n]))) < 1e-7
    end
end
