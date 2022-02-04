using WTP
using Test

@testset "Miller indices" begin
    sizes = (4, 4, 4)
    @test (1, 1, 1) == miller_to_standard(sizes, (0, 0, 0), (0, 0, 0))
    @test (3, 4, 2) == miller_to_standard(sizes, (-3, -1, 1), (1, 0, 0))

    @test (0, 0, 0) == standard_to_miller(sizes, (1, 1, 1), (0, 0, 0))
    @test (-3, -1, 1) == standard_to_miller(sizes, (3, 4, 2), (1, 0, 0))

    one_dim_index = three_to_one(miller_to_standard(sizes, (-3, -1, 1), (1, 0, 0))..., sizes)
    @test (-3, -1, 1) == standard_to_miller(sizes, one_to_three(one_dim_index, sizes), (1, 0, 0))
end