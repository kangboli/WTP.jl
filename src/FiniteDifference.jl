
export FiniteDifference,
    W90FiniteDifference,
    finite_difference,
    find_shells,
    compute_weights,
    spread,
    shells,
    weights,
    populate_integral_table!,
    center,
    neighbor_basis_integral,
    second_moment

abstract type FiniteDifference end

"""
    shells(scheme)

Shells of neighbors of the gamma point within a finite different scheme.
Each shell is a vector of kpoints.
"""
shells(scheme::FiniteDifference)::AbstractVector{Vector{KPoint}} =
    scheme.neighbor_shells

"""
    weights(scheme)

The weights corresponding to each shell within a scheme. The weights are ordered
from the inner-most to the outer-most shell.
"""
weights(scheme::FiniteDifference) = scheme.weights

neighbor_basis_integral(scheme::FiniteDifference) = scheme.neighbor_basis_integral
# transformed_integral(scheme::FiniteDifference) = scheme.transformed_integral

"""
    find_neighbors(k, scheme)

Find the set of relevant neighbors for a kpoint under a scheme.
"""
function find_neighbors(kpoint::KPoint, scheme::FiniteDifference)
    dk_list = vcat(shells(scheme)...)
    return (dk -> kpoint + dk).([dk_list; -dk_list])
end

"""
The finite difference scheme used in Wannier90.
"""
struct W90FiniteDifference <: FiniteDifference
    neighbor_shells::AbstractVector{Vector{KPoint}}
    weights::AbstractVector{Number}
    neighbor_basis_integral::NeighborIntegral
end

"""
The Brillouin zone on which the finite difference scheme is defined.
"""
grid(scheme::W90FiniteDifference) = grid(shells(scheme)[1][1])

function find_shells(u::Wannier, n_shell::Int)
    shells = SortedDict{Real,Vector{KPoint}}()
    d = collect(-2*n_shell:2*n_shell)
    for i in d, j in d, k in d
        k = KPoint(grid(u), [i, j, k], true)
        key = round(norm(cartesian(k)), digits = 5)
        haskey(shells, key) ? append!(shells[key], [k]) : shells[key] = [k]
    end

    # for k in grid(u)
    #     key = round(norm(cartesian(k)), digits = 5)
    #     haskey(shells, key) ? append!(shells[key], [k]) : shells[key] = [k]
    # end
    return collect(values(shells))[2:n_shell+1]
end

"""
Solve Aw = q
"""
function compute_weights(neighbor_shells::AbstractVector{Vector{KPoint}})

    indices = SortedDict(
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
        A[:, s] =
            [sum((b -> c(b, i) * c(b, j)).(neighbor_shells[s])) for (i, j) in keys(indices)]
    end

    q = [i == b ? 1 : 0 for (i, b) in keys(indices)]
    w = A \ q

    return isapprox(A * w, q, atol = 1e-13) ? w : nothing
end

function W90FiniteDifference(u::Wannier, n_shells = 1)
    neighbor_shells = find_shells(u, n_shells)
    weights = compute_weights(neighbor_shells)
    weights === nothing && return W90FiniteDifference(u, n_shells + 1)
    neighbor_integral = NeighborIntegral()

    scheme = W90FiniteDifference(neighbor_shells, weights, neighbor_integral)
    populate_integral_table!(scheme, u)
    return scheme
end

function center(M::NeighborIntegral, scheme::W90FiniteDifference, n::Int)
    function kpoint_contribution(k::KPoint)
        -sum(zip(weights(scheme), shells(scheme))) do (w, shell)
            sum(shell) do b
                w * cartesian(b) * imag(log(M[k, k+b][n, n]))
            end
        end
    end

    brillouin_zone = collect(grid(scheme))
    return sum(kpoint_contribution.(brillouin_zone)) / prod(size(brillouin_zone))
end

function second_moment(M::NeighborIntegral, scheme::W90FiniteDifference, n::Int)
    function kpoint_contribution(k::KPoint)
        sum(zip(weights(scheme), shells(scheme))) do (w, shell)
            sum(shell) do b
                w * (1 - abs2(M[k, k+b][n, n]) + imag(log(M[k, k+b][n, n]))^2)
            end
        end
    end
    brillouin_zone = collect(grid(scheme))
    return sum(kpoint_contribution.(brillouin_zone)) / prod(size(brillouin_zone))
end

function populate_integral_table!(scheme::W90FiniteDifference, wannier::Wannier)
    brillouin_zone = grid(wannier)
    M = neighbor_basis_integral(scheme)

    """
    Compute the mmn matrix between k1 and k2.
    """
    function mmn_matrix(k1::KPoint, k2::KPoint)::Matrix{ComplexFxx}
        return [braket(dagger(m), n) for m in wannier[k1], n in wannier[k2]]
    end

    # extended_brillouin_zone = union(Set{KPoint}([]), (k->find_neighbors(k, scheme)).(brillouin_zone)...)

    # Threads.@threads for kpoint in collect(brillouin_zone)
    @showprogress for kpoint in collect(brillouin_zone)
        # println(i, ", kpoint: ", coefficients(kpoint), ", n neighbor: ", length(neighbors))
        for neighbor in find_neighbors(kpoint, scheme)
            M[neighbor, kpoint] !== nothing && continue
            M[kpoint, neighbor] = mmn_matrix(kpoint, neighbor)
        end
    end
    for kpoint in brillouin_zone
        (m -> cache!(m, M)).(wannier[kpoint])
    end
    return M
end

