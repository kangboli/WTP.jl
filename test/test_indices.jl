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

## Test that miller_to_standard and standard_to_miller are inverse of each other.
# TODO: Test with offsets.

even_grid_miller = ((-3, 2), (-4, 3), (-2, 1))
odd_grid_miller = ((-3, 3), (-4, 4), (-2, 2))

for grid in [even_grid_miller, odd_grid_miller]
    ((x_1, x_2), (y_1, y_2), (z_1, z_2)) = grid
    sizes = x_2-x_1, y_2-y_1, z_2-z_1
    for i=x_1:x_2+1, j=y_1:y_2+1, k=z_1:z_2+1
        [i,j,k] == standard_to_miller(sizes, miller_to_standard(sizes, [i, j, k], [0, 0, 0]), [0,0,0])
    end
end