using WTP
using Test
using LinearAlgebra

@testset "Single WFC" begin
    wave_functions_filename = joinpath(test_1_dir, "si.save/wfc1.dat")
    wave_functions = WFC(wave_functions_filename)
    @test wave_functions.i_kpoint == 1
    @test wave_functions.k_coordinates == [0, 0, 0]
    @test !wave_functions.gamma_only
    @test wave_functions.n_planewaves == 537
    @test wave_functions.max_n_planewaves == 537
    @test wave_functions.n_band == 20
    @test size(wave_functions.miller) == (3, 537)
    @test wave_functions.evc === nothing

    load_evc!(wave_functions)
    @test size(wave_functions.evc) == (537, 20)
end

@testset "Contructing the Brillouin Zone" begin
    wave_functions_list = wave_functions_from_directory(joinpath(test_1_dir, "si.save"))
    coordinates_kpoints = [wave_functions.k_coordinates for wave_functions in wave_functions_list]
    brillouin_zone = brillouin_zone_from_k_coordinates(coordinates_kpoints, wave_function_basis(wave_functions_list[1]))
    @test size(brillouin_zone) == (4, 4, 4)
    wannier = init_orbital_set(brillouin_zone)
    sizes = Tuple(maximum((wave_function) -> estimate_sizes(wave_function, i, 2), wave_functions_list) for i in 1:3)
    @test sizes == (24, 24, 24)
end

@testset "Orbital" begin
    wave_functions_list = wave_functions_from_directory(joinpath(test_1_dir, "si.save"))
    ũ = orbital_set_from_save(wave_functions_list);
    gamma_point = grid(ũ)[0, 0, 0]
    u1gamma = ũ[gamma_point][1]
    @test size(grid(u1gamma)) == (24, 24, 24)
    @test ket(u1gamma) == true
    # Inner product should give the unit norm.
    @test isapprox(abs(braket(dagger(u1gamma), u1gamma)), 1.0, atol=1e-7)
    @test !ket(dagger(u1gamma))
end