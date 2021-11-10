@test "spread" begin
    wannier = wannier_from_save(joinpath(experiment_dir, "si.save"));
    g = grid(wannier)
    integrals = construct_integral_table(wannier)
    # integrals[g[0,0,0]][g[0,0,2]]
    integrals[g[0,0,0]][g[0,0,1]]
    integrals[g[0,0,0], g[1,0,0]][5, 5]

    integrals[g[0,0,0], g[0,0,1]][5, 5]

    integrals[g[0,0,0], g[0,-1,0]][5, 5]
    integrals[g[0,0,0], g[0,-1,0]][5, 5]

    integrals[g[0,0,0], g[1,0,0]][2, 2]
    integrals[g[0,0,0], g[1,0,0]][3, 3]
    integrals[g[-1,0,0], g[0,-1,0]]
    integrals[g[3,2,-1], g[2,3,-1]]
    # integrals[g[0,0,0], g[-2,0,0]]

    # braket(dagger(wannier[g[0,0,0]], 1), wannier[g[0,0,2]])


    @time gradient = finite_difference(wannier, 1, g[0,0,0])
    n = 2
    graident = finite_difference(wannier, n, g[0,0,0])
    [braket(dagger(wannier[n, g[0,0,0]]), gradient[i]) for i = 1:3]
    
    braket(dagger(wannier[1, g[0,0,0]]), gradient[2])
    braket(dagger(wannier[1, g[0,0,0]]), gradient[3])

braket(dagger(wannier[1, g[0,0,0]]), gradient[1])
    s, c = spread(wannier, 1)
    s, c = spread(wannier, 2)
    s, c = spread(wannier, 5)
    s, c = spread(wannier, 7)
    s - norm(c)^2
end


wannier = wannier_from_unk_dir(joinpath(e01_dir, "unk"),
                               joinpath(e01_dir, "si.save"))
scheme = W90FiniteDifference(wannier, 1)
integrals = construct_integral_table(wannier, scheme)

s, c = spread(wannier, 1, scheme)
s, c = spread(wannier, 2, scheme)
s, c = spread(wannier, 13, scheme)
s, c = spread(wannier, 4, scheme)

# for k in grid(wannier)
#     println(length(wannier[k]))
# end

# s, c = spread(wannier, 1)
# s, c = spread(wannier, 6)
# integrals[g[0,0,0], g[1,0,0]][1, 1]
# integrals[g[0,0,0], g[0,0,1]][1, 1]
# b_1, b_2, b_3 = basis(grid(g[1,0,0]))
# collect(grid(g[1,0,0]))
# B = vector3_to_matrix(basis(g))
# norm(coefficients(g[1,0,0]))
# 3 / (8*norm(basis_transform(coefficients(g[1,0,0]), basis(g[1,0,0]), CARTESIAN_BASIS))^2)
# inv(B)


# neighbor_shells = find_shells(wannier, 1)
# weights = WTP.compute_weights(neighbor_shells)

# grad = finite_difference(wannier, 1, g[1,0, 0], neighbor_shells, weights)
