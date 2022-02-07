using WTP
using Test

@testset "Miller indices" begin

    sizes = (3, 5, 6)
    A = rand(Float64, sizes)

    ## Test conversion from one index to three.
    for i = 1:prod(sizes)
        @test A[i] == A[one_to_three(i, sizes)...]
    end

    ## Test conversion from three indices to one.
    for i=1:3
        for j=1:5
            for k=1:6
                @test A[i, j, k] == A[three_to_one(i, j, k, sizes)]
            end
        end
    end

    ## Test reshape.
    for i = 1:prod(sizes)
        @test reshape(A, prod(sizes))[i] == A[one_to_three(i, sizes)...]
    end


    sizes = (4, 4, 4)
    @test (1, 1, 1) == miller_to_standard(sizes, (0, 0, 0), (0, 0, 0))
    @test (3, 4, 2) == miller_to_standard(sizes, (-3, -1, 1), (1, 0, 0))
    @test (0, 0, 0) == standard_to_miller(sizes, (1, 1, 1), (0, 0, 0))
    @test (-3, -1, 1) == standard_to_miller(sizes, (3, 4, 2), (1, 0, 0))

    one_dim_index = three_to_one(miller_to_standard(sizes, (-3, -1, 1), (1, 0, 0))..., sizes)
    @test (-3, -1, 1) == standard_to_miller(sizes, one_to_three(one_dim_index, sizes), (1, 0, 0))
    ## Test that miller_to_standard and standard_to_miller are inverse of each other.
    # TODO: Test with offsets.

    even_grid_miller = ((-3, 2), (-4, 3), (-2, 1))
    odd_grid_miller = ((-3, 3), (-4, 4), (-2, 2))

    for grid in [even_grid_miller, odd_grid_miller]
        ((x_1, x_2), (y_1, y_2), (z_1, z_2)) = grid
        sizes = (x_2 - x_1+1, y_2 - y_1+1, z_2 - z_1+1)
        for i = x_1:x_2, j = y_1:y_2, k = z_1:z_2
            standard_indices = miller_to_standard(sizes, (i, j, k), (0, 0, 0))
            @test (i, j, k) == standard_to_miller(sizes, standard_indices, (0, 0, 0))
        end
    end
end