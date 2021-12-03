
export FiniteDifference,
    W90FiniteDifference,
    W90FiniteDifference3D,
    finite_difference,
    find_shells,
    compute_weights,
    spread,
    shells,
    weights,
    populate_integral_table!,
    center,
    neighbor_basis_integral,
    second_moment,
    BranchStable

abstract type FiniteDifference end

"""
    shells(scheme)

Shells of neighbors of the gamma point within a finite different scheme.
Each shell is a vector of kpoints.
"""
shells(scheme::FiniteDifference)::AbstractVector{Vector{<:KPoint}} =
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
    # TODO: The negative part may not be necessary.
    return (dk -> kpoint + dk).([dk_list; -dk_list])
end

abstract type W90FiniteDifference <: FiniteDifference end

"""
The finite difference scheme used in Wannier90.
"""
struct W90FiniteDifference3D <: W90FiniteDifference
    neighbor_shells::AbstractVector{Vector{<:KPoint}}
    weights::AbstractVector{Number}
    neighbor_basis_integral::NeighborIntegral
end

"""
The Brillouin zone on which the finite difference scheme is defined.
"""
grid(scheme::FiniteDifference) = grid(shells(scheme)[1][1])

function find_shells(grid::Grid, n_shell::Int)
    shells = SortedDict{Real,Vector{KPoint}}()
    d = collect(-2*n_shell:2*n_shell)
    for i in d, j in d, k in d
        k = grid_vector_constructor(grid, [i, j, k])
        key = round(norm(cartesian(k)), digits = 5)
        haskey(shells, key) ? append!(shells[key], [k]) : shells[key] = [k]
    end

    return collect(values(shells))[2:n_shell+1]
end

"""
Solve Aw = q
"""
function compute_weights(neighbor_shells::Vector{Vector{T}}) where T <: AbstractGridVector

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

function W90FiniteDifference3D(u::Wannier{UnkBasisOrbital{ReciprocalLattice3D}}, n_shells = 1)
    neighbor_shells = find_shells(grid(u), n_shells)
    weights = compute_weights(neighbor_shells)
    weights === nothing && return W90FiniteDifference3D(u, n_shells + 1)
    neighbor_integral = NeighborIntegral()

    scheme = W90FiniteDifference3D(neighbor_shells, weights, neighbor_integral)
    populate_integral_table!(scheme, u)
    return scheme
end

abstract type BranchStable end

function center(M::NeighborIntegral, scheme::W90FiniteDifference, n::Int)
    function kpoint_contribution(k::KPoint)
        -sum(zip(weights(scheme), shells(scheme))) do (w, shell)
            sum(shell) do b
                w * cartesian(b) * angle(M[k, k+b][n, n])
            end
        end
    end

    brillouin_zone = collect(grid(scheme))
    return sum(kpoint_contribution.(brillouin_zone)) / prod(size(brillouin_zone))
end

function center(M::NeighborIntegral, scheme::W90FiniteDifference, n::Int, ::Type{BranchStable})
    function kpoint_contribution(k::KPoint)
        -sum(zip(weights(scheme), shells(scheme))) do (w, shell)
            sum(unique(k->Set([k, -k]), shell)) do b
                ϕ⁺ = M[k, k+b][n, n] |> angle
                ϕ⁻ = M[k, k-b][n, n] |> angle
                branch = (sign(ϕ⁺) == sign(ϕ⁻) ? -1 : 1)
                w * cartesian(b) * ϕ⁺ + w * cartesian(-b) * branch * ϕ⁻
            end
        end
    end

    brillouin_zone = scheme |> grid |> collect
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

function populate_integral_table!(scheme::FiniteDifference, u::Wannier)
    brillouin_zone = grid(u)
    M = neighbor_basis_integral(scheme)

    """
    Compute the mmn matrix between k1 and k2.
    """
    function mmn_matrix(k_1::T, k_2::T) where T <: AbstractGridVector{<:BrillouinZone}
        return [braket(dagger(m), n) for m in u[k_1], n in u[k_2]]
    end

    # Threads.@threads for kpoint in collect(brillouin_zone)
    @showprogress for kpoint in collect(brillouin_zone)
        for neighbor in find_neighbors(kpoint, scheme)
            M[neighbor, kpoint] !== nothing && continue
            M[kpoint, neighbor] = mmn_matrix(kpoint, neighbor)
        end
    end
    for kpoint in brillouin_zone
        (m -> cache!(m, M)).(u[kpoint])
    end
    return M
end


# extended_brillouin_zone = union(Set{KPoint}([]), (k->find_neighbors(k, scheme)).(brillouin_zone)...)