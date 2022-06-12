using WTP
using Test
using LinearAlgebra

@testset "Convolution" begin

    ## Load the wave function.
    wave_functions_list = wave_functions_from_directory(joinpath(test_6_dir, "si.save"))
    ũ = orbital_set_from_save(wave_functions_list)
    brillouin_zone = grid(ũ)

    ## Get the map of kpoints.
    k_map, _ = i_kpoint_map(wave_functions_list)


    ## Load the gauge transform amn and verify it.
    amn = AMN(joinpath(test_6_dir, "output/pw2wan/si.amn"))
    U = Gauge(grid(ũ), amn, k_map)
    ũ = set_gauge(ũ, U)
    u = ifft(ũ)

    supercell_ = expand(transform_grid(orbital_grid(ũ)), [size(brillouin_zone)...])
    r̃2 = fft(map(r -> norm(r)^2, supercell_), false)
    wannier_functions = commit_gauge(ũ)(:)
    densities = abs2.(ifft.(wannier_functions))

    # Verify that the densities calculated efficiently agree with that calculated more patiently.
    phase = phase_factors(u)
    densities_slow = abs2.([u(i, phase) for i = 1:4])

    for i = 1:4
        @test isapprox(norm(elements(densities_slow[i] - densities[i])), 0, atol=1e-7)
    end
    
    ρ̃_list = (ρ -> fft(ρ, false)).(densities)

    # Construct the scheme and apply the transform.
    # scheme = CosScheme3D(ũ)
    # M = gauge_transform(neighbor_basis_integral(scheme), U)

    # Check that the center and the spread are not too far away from W90.
    for i = 1:4
        c, σ = center_spread(reciprocal_densities[i], r̃2)
        # println(c)
        # println(σ)
        for j = 1:3
            @test isapprox(abs(c[j]), 1.3, atol=0.5)
        end
        @test isapprox(σ, 7, atol=0.5)
        # @test isapprox(norm(c - center(M, scheme, i, BranchStable)), 0, atol=2e-2)
        # @test isapprox(σ - spread(M, scheme, i), 0, atol=2)
    end
end