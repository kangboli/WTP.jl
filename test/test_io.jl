using WTP
using Test
using LinearAlgebra

@testset "AMN and MMN" begin
    
    # Read 
    amn = AMN(joinpath(test_2_dir, "output/pw2wan/si.amn"))
    # Write
    WTP.dump(amn, joinpath(test_2_dir, "si_dumped.amn"))
    # Read again and make sure that they are the same.
    amn_2 = AMN(joinpath(test_2_dir, "si_dumped.amn"))
    @test amn_2.gauge == amn.gauge
    @test amn.n_band == amn_2.n_band
    @test amn.n_kpoint == amn_2.n_kpoint
    @test amn.n_wannier == amn_2.n_wannier

    wave_functions_list = wave_functions_from_directory(joinpath(test_2_dir, "si.save"))
    k_map, brillouin_zone = i_kpoint_map(wave_functions_list)
    U = Gauge(brillouin_zone, amn, k_map, false)
    amn_3 = AMN(U, k_map)

    @test amn_3.gauge == amn.gauge
    @test amn_3.n_band == amn.n_band
    @test amn_3.n_kpoint == amn.n_kpoint
    @test amn_3.n_wannier == amn.n_wannier
end