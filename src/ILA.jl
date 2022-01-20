
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
    second_moment,
    spread,
    neighbor_basis_integral,
    BranchStable,
    BranchNaive,
    ila_gradient,
    optimize

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

function gauge_transform(scheme::FiniteDifference, U::Gauge)
    @set scheme.neighbor_basis_integral = gauge_transform(neighbor_basis_integral(scheme), U)
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
abstract type BranchNaive end

center(scheme::W90FiniteDifference, n::Integer) =center(scheme, n, BranchNaive) 
center(scheme::W90FiniteDifference, n::Integer, ::Type{T}) where T = center(neighbor_basis_integral(scheme), scheme, n, T)
center(M::NeighborIntegral, scheme::W90FiniteDifference, n::Integer) = center(M, scheme, n, BranchNaive)

function center(M::NeighborIntegral, scheme::W90FiniteDifference, n::Int, ::Type{BranchNaive})
    function kpoint_contribution(k::KPoint)
        -sum(zip(weights(scheme), shells(scheme))) do (w, shell)
            sum(shell) do b
                w * cartesian(b) * angle(M[k, k+b][n, n])
            end
        end
    end

    brillouin_zone = collect(grid(scheme))
    return sum(kpoint_contribution.(brillouin_zone)) / length(brillouin_zone)
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

second_moment(scheme::W90FiniteDifference, n::Int) = second_moment(neighbor_basis_integral(scheme), scheme, n)

function second_moment(M::NeighborIntegral, scheme::W90FiniteDifference, n::Int)
    function kpoint_contribution(k::KPoint)
        sum(zip(weights(scheme), shells(scheme))) do (w, shell)
            sum(shell) do b
                w * (1 - abs2(M[k, k+b][n, n]) + angle(M[k, k+b][n, n])^2)
            end
        end
    end
    brillouin_zone = collect(grid(scheme))
    return sum(kpoint_contribution.(brillouin_zone)) / prod(size(brillouin_zone))
end

spread(scheme::W90FiniteDifference, n::Integer) = spread(scheme, n, BranchNaive)
spread(scheme::W90FiniteDifference, n::Integer, ::Type{T}) where T = spread(neighbor_basis_integral(scheme), scheme, n, T)
spread(M::NeighborIntegral, scheme::W90FiniteDifference, n::Integer) = spread(M, scheme, n, BranchNaive) 
spread(M::NeighborIntegral, scheme::W90FiniteDifference, n::Integer, ::Type{T}) where T = second_moment(M, scheme, n) - 
    norm(center(M, scheme, n, T))^2
    


function populate_integral_table!(scheme::FiniteDifference, u::Wannier)
    brillouin_zone = grid(u)
    M = neighbor_basis_integral(scheme)

    """
    Compute the mmn matrix between k1 and k2.
    """
    function mmn_matrix(k_1::T, k_2::T) where T <: AbstractGridVector{<:BrillouinZone}

        U_adjoint = vcat(adjoint.(vectorize.(u[k_1]))...)
        V = hcat(map(vectorize, u[k_2])...)
        return U_adjoint * V
    end

    @showprogress for k in collect(brillouin_zone)
        for neighbor in find_neighbors(k, scheme)
            M[neighbor, k] !== nothing && continue
            # M[k, neighbor] = adjoint(U[k]) * mmn_matrix(k, neighbor) * U[neighbor]
            M[k, neighbor] = mmn_matrix(k, neighbor)
        end
    end
    for k in brillouin_zone
        (m -> cache!(m, M)).(u[k])
    end
    return M
end

function n_band(M::NeighborIntegral)
    first_matrix = collect(values(integrals(M)))[1]
    return size(first_matrix, 1)
end

function ila_gradient(M::NeighborIntegral, scheme::W90FiniteDifference, brillouin_zone::B) where B <: BrillouinZone
    N = n_band(M)
    c = (n->center(M, scheme, n)).(1:N)
    
    G = k-> sum(zip(weights(scheme), shells(scheme))) do (w, shell)
        sum(shell) do b
            A = M[k, k+b]
            q = [angle(A[n, n]) + cartesian(b)' * c[n] for n=1:N]
            R = hcat([A[:, n] * A[n, n]' for n=1:N]...)
            T = hcat([(A[:, n] / A[n, n]) * q[n] for n=1:N]...)
            4w * ((R-R')/2 - (T+T')/2im) / length(brillouin_zone)
        end
    end

    return map(G, brillouin_zone)
end


total_spread(M, scheme) = let N = n_band(M)
    sum(1:N) do n
        spread(M, scheme, n)
    end
end

function optimize(U_0::Gauge, scheme::W90FiniteDifference, brillouin_zone::B, ϵ=1e-7) where B <: BrillouinZone
    α = 2
    M_0 = neighbor_basis_integral(scheme)
    N = n_band(M_0)
    U = U_0
    M = gauge_transform(M_0, U)
    Ω = total_spread(M, scheme)

    function line_search(U, Ω, ∇Ω, ∇Ω²)
        while true
            α /= 2
            Q = Ω - 0.5α * ∇Ω²
            new_elements = elements(map(k->U[k] * cis(-Hermitian(1im * α * ∇Ω[k])), brillouin_zone))
            V = Gauge{B}(brillouin_zone, new_elements)
            M_V = gauge_transform(M_0, V)
            Ω_V = total_spread(M_V, scheme)
            println("Q: $(Q)"); println("Ω: $(Ω)"); println("Ω_V: $(Ω_V)"); println("α: $(α)"); flush(stdout)
            Ω_V > Q && continue
            α *= 2
            return V, M_V, Ω_V
        end
    end

    while true
        ∇Ω = Gauge{B}(brillouin_zone, elements(ila_gradient(M, scheme, brillouin_zone))) 
        ∇Ω² = sum(brillouin_zone) do k
            v = reshape(∇Ω[k], N^2)
            abs(v' * v)
        end
        println("∇Ω²: $(∇Ω²)")
        √∇Ω² <= ϵ && break
        U, M, Ω = line_search(U, Ω, ∇Ω, ∇Ω²)
    end

    return U
end