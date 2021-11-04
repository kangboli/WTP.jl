export FiniteDifference, W90FiniteDifference, finite_difference, find_shells,
compute_weights, spread, neighbor_shells, weights

abstract type FiniteDifference end

neighbor_shells(scheme::FiniteDifference) = scheme.neighbor_shells
weights(scheme::FiniteDifference) = scheme.weights

function find_neighbors(kpoint::KPoint, scheme::FiniteDifference)
    dk_list = vcat(neighbor_shells(scheme)...)
    return (dk -> kpoint + dk).([dk_list; -dk_list])
end

struct W90FiniteDifference <: FiniteDifference
    neighbor_shells::AbstractVector{Vector{KPoint}}
    weights::AbstractVector{Number}
end


function find_shells(u::Wannier, n_shell::Int)
    shells = OrderedDict{Float64,Vector{KPoint}}()
    
    for k in grid(u)
        key = round(norm(cartesian(k)), digits = 7)
        haskey(shells, key) ? append!(shells[key], [k]) : shells[key] = [k]
    end
    return collect(values(shells))[1:n_shell]
end

function compute_weights(neighbor_shells::AbstractVector{Vector{KPoint}})

    indices = OrderedDict(
        (1, 1) => 1,
        (1, 2) => 2,
        (1, 3) => 3,
        (2, 2) => 4,
        (2, 3) => 5,
        (3, 3) => 6,
    )

    A = zeros((6, length(neighbor_shells)))
    c(b, i) = cartesian(b)[i]
    for s = 1:length(neighbor_shells)
        A[:, s] = [sum((b -> c(b, i) * c(b, j)).(neighbor_shells[s])) for (i, j) in keys(indices)]
    end

    q = [i == b ? 1 : 0 for (i, b) in keys(indices)]
    w = A \ q

    return isapprox(A * w, q, atol = 1e-13) ? w : nothing
end

function W90FiniteDifference(u::Wannier, n_shells = 1)
    neighbor_shells = find_shells(u, n_shells)
    weights = compute_weights(neighbor_shells)
    return weights === nothing ? W90FiniteDifference(u, n_shells + 1) :
           W90FiniteDifference(neighbor_shells, weights)
end


"""
"""
function finite_difference(wannier::Wannier, n::Int, k::KPoint, scheme::W90FiniteDifference)
    f(k) = wannier[n, k]
    grad = zeros(UnkOrbital, 3)
    for (w, shell) in zip(weights(scheme), neighbor_shells(scheme))
        grad += sum((b -> w * cartesian(b) * (f(k + b) - f(k))).(shell))
    end

    return grad
end


function spread(u::Wannier, n::Int64, scheme::FiniteDifference)

    center = zeros(ComplexF64, 3)
    r2::ComplexF64 = 0
    N = prod(size(grid(u)))
    for k in grid(u)
        # println("adding $(coefficients(k))")
        # gradient = finite_difference(u, n, k, scheme)
        # center += [(1im / N) * braket(dagger(u[n, k]), g) for g in gradient]
        for (w, shell) in zip(weights(scheme), neighbor_shells(scheme))
            for b in shell
                center += - 1 / N * w * cartesian(b) * imag(log(braket(dagger(u[n, k]), u[n, k+b])))
            end
        end

        # println(center)
        # r2 += (1 / N) * sum([braket(dagger(gradient[i]), gradient[i]) for i = 1:3])
    end
    return r2, center
end

function construct_integral_table(wannier::Wannier, scheme::FiniteDifference)
    brillouin_zone = grid(wannier)
    neighbor_set = NeighborIntegral()

    """
    Compute the mmn matrix between k1 and k2.
    """
    function mmn_matrix(k1::KPoint, k2::KPoint)::Matrix{ComplexF64}
        return [braket(dagger(m), n) for m in wannier[k1], n in wannier[k2]]
    end

    extended_brillouin_zone = union(Set{KPoint}([]), (k->find_neighbors(k, scheme)).(brillouin_zone)...)
    i = Threads.Atomic{Int}(0)

    # Threads.@threads for kpoint in collect(brillouin_zone)
    for kpoint in extended_brillouin_zone
        neighbors = Set(find_neighbors(kpoint, scheme))
        println(i, ", kpoint: ", coefficients(kpoint), ", n neighbor: ", length(neighbors))
        for neighbor in filter(x -> x in extended_brillouin_zone, find_neighbors(kpoint, scheme))
            neighbor_set[neighbor, kpoint] !== nothing && continue
            neighbor_set[kpoint, neighbor] = mmn_matrix(kpoint, neighbor)
        end
        Threads.atomic_add!(i, 1)
    end
    for kpoint in brillouin_zone
        (m -> cache!(m, neighbor_set)).(wannier[kpoint])
    end
    return neighbor_set
end


# function find_neighbors(kpoint::KPoint)
#     b1 = grid(kpoint)[1, 0, 0]
#     b2 = grid(kpoint)[0, 1, 0]
#     b3 = grid(kpoint)[0, 0, 1]
#     dk_list = [b1, b2, b3, b1 - b2, b2 - b3, b3 - b1]
#     # Second order list.
#     # dk_list = [ 
#     # b1+b1, b1+b2, b1+b3,
#     # b2+b1, b2+b2, b2+b3,
#     # b3+b1, b3+b2, b3+b3,
#     # b1-b2, b1-b3,
#     # b2-b1, b2-b3,
#     # b3-b1, b3-b2,
#     # ]
#     return (dk -> kpoint + dk).([dk_list; -dk_list])
# end
