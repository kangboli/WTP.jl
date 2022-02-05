
using WTP
using Test
using LinearAlgebra

@testset "Truncated Convolution" begin
    ## Load the wave function.

    wave_functions_list = wave_functions_from_directory(joinpath(test_6_dir, "si.save"))
    ũ = wannier_from_save(wave_functions_list)
    brillouin_zone = grid(ũ)

    ## Get the map of kpoints.

    k_map, _ = i_kpoint_map(wave_functions_list)

    ## Load the gauge transform amn and verify it.
    amn = AMN(joinpath(test_6_dir, "output/pw2wan/si.amn"))
    U = Gauge(grid(ũ), amn, k_map)

    ## Construct the scheme and compare the gradients.
    scheme = W90FiniteDifference3D(ũ)
    G_1 = gauge_gradient(U, scheme, brillouin_zone, TruncatedConvolution)
    G_2 = gauge_gradient(U, scheme, brillouin_zone, BranchNaive)

    @test sum(k -> norm(G_1[k] - G_2[k]), brillouin_zone) < 1

    ## Apply the gauge transform on the basis integrals to get the neighbor integrals.

    M = gauge_transform(neighbor_basis_integral(scheme), U)

    ## Get the center from the neighbor integrals and the scheme. Test it against the output 
    # from Wannier90.

    for i = 1:4
        c_1 = center(M, scheme, i, TruncatedConvolution)
        s_1 = spread(M, scheme, i, TruncatedConvolution)
        c_2 = center(M, scheme, i, BranchNaive)
        s_2 = spread(M, scheme, i, BranchNaive)
        @test isapprox(norm(c_1 - c_2), 0, atol = 1e-7)
        @test isapprox(s_1 - s_2, 0, atol = 0.5)
    end


end