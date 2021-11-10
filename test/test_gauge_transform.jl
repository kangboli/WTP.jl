wave_functions_list = wave_functions_from_directory(joinpath(test_2_dir, "si.save"))
wannier = wannier_from_save(wave_functions_list);
g = grid(wannier)
k_map, _ = i_kpoint_map(wave_functions_list)

# Load the gauge transform amn and verify it.
amn = AMN(joinpath(test_2_dir, "output/pw2wan/si.amn"))
U = Gauge(grid(wannier), amn, k_map)

# These numbers are taken from the amn file.
isapprox(norm(U[g[0, 0, 0]][1, 1] - (0.220397387839 - 0.046950436272im)), atol=1e-7)
isapprox(norm(U[g[0, 0, 0]][2, 1] - (0.207291128273 - 0.299289327186im)), atol=1e-7)


# Load the neighbor integrals mmn. 
mmn = MMN(joinpath(test_2_dir, "output/pw2wan/si.mmn"))
neighbor_integral = NeighborIntegral(mmn, k_map)

# Perform the gauge transform and verify that the result is correct.
neighbor_integral = gauge_transform(neighbor_integral, U)

d = diag(neighbor_integral[g[0, 0, 0], g[0, 1, 0]])
# These numbers are obtained by printing the m_matrix from Wannier90.
@test isapprox(real(d[1]), 0.51130606080639185, atol = 1e-6)
@test isapprox(imag(d[1]), 0.54468481207629793, atol = 1e-6)
@test isapprox(real(d[2]), 0.68760248170789728, atol = 1e-6)
@test isapprox(imag(d[2]), -0.24726067511852651, atol = 1e-6)
@test isapprox(real(d[3]), 0.95497390215381983, atol = 1e-6)
@test isapprox(imag(d[3]), 2.14367698922523602E-002, atol = 1e-6)


scheme = W90FiniteDifference(wannier)

n, N = 3, length(collect(grid(wannier)))
center = zeros(ComplexFxx, 3)

for k in grid(wannier)
    for (w, shell) in zip(weights(scheme), neighbor_shells(scheme))
        for b in shell
            # neighbor_integral[k, k+b]
            center += -1 / N * w * cartesian(b) * imag(log(neighbor_integral[k, k+b][n, n]))
        end
    end
end

center


for ((k_1, k_2), integral) in integrals(neighbor_integral)
    first = findfirst(
        x -> isapprox(abs(x - 0.5113060608063918 + 0.54468481207629793im), 0, atol = 1e-5),
        integral,
    )
    first !== nothing && println(k1, k2)
end

# @test isapprox(real(d[4]), 0.94609710104871247, atol = 1e-6)
# @test isapprox(imag(d[4]), 2.24223892764866262E-002, atol = 1e-6)
# @test isapprox(real(d[5]), 0.86381202162192472, atol = 1e-6)
# @test isapprox(imag(d[5]), -3.61638924620649033E-002, atol = 1e-6)
# @test isapprox(real(d[6]), 0.82762311368765762, atol = 1e-6)
# @test isapprox(imag(d[6]), 0.35711407338455936, atol = 1e-6)
# @test isapprox(real(d[7]), 0.85524734045560435, atol = 1e-6)
# @test isapprox(imag(d[7]), 0.16399655196515897, atol = 1e-6)
# @test isapprox(real(d[8]), 0.88334677577137455, atol = 1e-6)
# @test isapprox(imag(d[8]), 0.18720486670427361, atol = 1e-6)
# @test isapprox(real(d[9]), 0.89201457827146113, atol = 1e-6)
# @test isapprox(imag(d[9]), 0.14457887209362333, atol = 1e-6)
# @test isapprox(real(d[10]), 0.71436047088676602, atol = 1e-6)
# @test isapprox(imag(d[10]), 0.54146577430624077, atol = 1e-6)
# @test isapprox(real(d[11]), 0.79191318540753108, atol = 1e-6)
# @test isapprox(imag(d[11]), 0.42138481733203798, atol = 1e-6)
# @test isapprox(real(d[12]), 0.88358878403767938, atol = 1e-6)
# @test isapprox(imag(d[12]), -0.21707809545708487, atol = 1e-6)
# @test isapprox(real(d[13]), 0.92716031970141899, atol = 1e-6)
# @test isapprox(imag(d[13]), 0.18861301625908158, atol = 1e-6)
# @test isapprox(real(d[14]), 0.72825992993729227, atol = 1e-6)
# @test isapprox(imag(d[14]), 0.49060227150616376, atol = 1e-6)
# @test isapprox(real(d[15]), 0.84928142191700173, atol = 1e-6)
# @test isapprox(imag(d[15]), -6.85431560648076810E-002, atol = 1e-6)
# @test isapprox(real(d[16]), 0.85486926865436108, atol = 1e-6)
# @test isapprox(imag(d[16]), 0.16796782496905474, atol = 1e-6)
# @test isapprox(real(d[17]), 0.77852576478582236, atol = 1e-6)
# @test isapprox(imag(d[17]), -0.56149595555280707, atol = 1e-6)
# @test isapprox(real(d[18]), 0.82733894196890878, atol = 1e-6)
# @test isapprox(imag(d[18]), -0.49778873189516665, atol = 1e-6)
# @test isapprox(real(d[19]), 0.71263308820421700, atol = 1e-6)
# @test isapprox(imag(d[19]), -0.61103172441259368, atol = 1e-6)
# @test isapprox(real(d[20]), 0.90618893074321860, atol = 1e-6)
# @test isapprox(imag(d[20]), -0.40473376765639657, atol = 1e-6)
