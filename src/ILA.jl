
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
    TruncatedConvolution,
    gauge_gradient,
    ILAOptimizer

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
function compute_weights(neighbor_shells::Vector{Vector{T}}) where {T<:AbstractGridVector}

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
abstract type TruncatedConvolution end

center(scheme::W90FiniteDifference, n::Integer) = center(scheme, n, BranchNaive)
center(scheme::W90FiniteDifference, n::Integer, ::Type{T}) where {T} = center(neighbor_basis_integral(scheme), scheme, n, T)
center(M::NeighborIntegral, scheme::W90FiniteDifference, n::Integer) = center(M, scheme, n, BranchNaive)

function center(M::NeighborIntegral, scheme::W90FiniteDifference, n::Int, ::Type{TruncatedConvolution})
    brillouin_zone = collect(grid(scheme))

    phase(b::KPoint) = angle(1 / length(brillouin_zone) * sum(k -> M[k, k+b][n, n], brillouin_zone))

    return -sum(zip(weights(scheme), shells(scheme))) do (w, shell)
        sum(unique(k -> Set([k, -k]), shell)) do b
            ϕ⁺, ϕ⁻ = phase(b), phase(-b)
            branch = (sign(ϕ⁺) == sign(ϕ⁻) ? -1 : 1)
            w * cartesian(b) * ϕ⁺ + w * cartesian(-b) * branch * ϕ⁻
        end
    end
end

function center(M::NeighborIntegral, scheme::W90FiniteDifference, n::Int, ::Type{BranchNaive})
    kpoint_contribution(k::KPoint) = -sum(zip(weights(scheme), shells(scheme))) do (w, shell)
        sum(b->w * cartesian(b) * angle(M[k, k+b][n, n]), shell) 
    end

    brillouin_zone = collect(grid(scheme))
    return sum(kpoint_contribution.(brillouin_zone)) / length(brillouin_zone)
end

function center(M::NeighborIntegral, scheme::W90FiniteDifference, n::Int, ::Type{BranchStable})
    kpoint_contribution(k::KPoint) = -sum(zip(weights(scheme), shells(scheme))) do (w, shell)
        sum(unique(k -> Set([k, -k]), shell)) do b
            ϕ⁺ = M[k, k+b][n, n] |> angle
            ϕ⁻ = M[k, k-b][n, n] |> angle
            branch = (sign(ϕ⁺) == sign(ϕ⁻) ? -1 : 1)
            w * cartesian(b) * ϕ⁺ + w * cartesian(-b) * branch * ϕ⁻
        end
    end

    brillouin_zone = scheme |> grid |> collect
    return sum(kpoint_contribution.(brillouin_zone)) / prod(size(brillouin_zone))

end

second_moment(scheme::W90FiniteDifference, n::Int) = second_moment(neighbor_basis_integral(scheme), scheme, n)

function second_moment(M::NeighborIntegral, scheme::W90FiniteDifference, n::Int)
    kpoint_contribution(k::KPoint) = sum(zip(weights(scheme), shells(scheme))) do (w, shell)
        sum(b->w * (1 - abs2(M[k, k+b][n, n]) + angle(M[k, k+b][n, n])^2), shell)
    end
    brillouin_zone = collect(grid(scheme))
    return sum(kpoint_contribution.(brillouin_zone)) / prod(size(brillouin_zone))
end

spread(scheme::W90FiniteDifference, n::Integer) = spread(scheme, n, BranchNaive)
spread(scheme::W90FiniteDifference, n::Integer, ::Type{T}) where {T} = spread(neighbor_basis_integral(scheme), scheme, n, T)
spread(M::NeighborIntegral, scheme::W90FiniteDifference, n::Integer) = spread(M, scheme, n, BranchNaive)
spread(M::NeighborIntegral, scheme::W90FiniteDifference, n::Integer, ::Type{T}) where {T} = second_moment(M, scheme, n) -
                                                                                            norm(center(M, scheme, n, T))^2

function spread(M::NeighborIntegral, scheme::W90FiniteDifference, n::Integer, ::Type{TruncatedConvolution})
    brillouin_zone = collect(grid(scheme))
    ρ̃(b) = sum(k -> M[k, k+b][n, n], brillouin_zone) / length(brillouin_zone) 

    sum(zip(weights(scheme), shells(scheme))) do (w, shell)
        sum(b -> 2w * (1 - abs(ρ̃(b))), shell)
    end
end


function populate_integral_table!(scheme::FiniteDifference, u::Wannier)
    brillouin_zone = grid(u)
    M = neighbor_basis_integral(scheme)

    """
    Compute the mmn matrix between k1 and k2.
    """
    function mmn_matrix(k_1::T, k_2::T) where {T<:AbstractGridVector{<:BrillouinZone}}
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
    return M
end

# for k in brillouin_zone
#     (m -> cache!(m, M)).(u[k])
# end

function n_band(M::NeighborIntegral)
    first_matrix = collect(values(integrals(M)))[1]
    return size(first_matrix, 1)
end

gauge_gradient(U::Gauge, scheme::W90FiniteDifference, brillouin_zone::BrillouinZone) =
    gauge_gradient(U::Gauge, scheme::W90FiniteDifference, brillouin_zone::BrillouinZone, BranchNaive)


"""
Gauge gradient for the finite difference.

"""
function gauge_gradient(U::Gauge, scheme::W90FiniteDifference, brillouin_zone::B, ::Type{BranchNaive}) where {B<:BrillouinZone}
    M = gauge_transform(neighbor_basis_integral(scheme), U)
    N = n_band(M)
    c = (n -> center(M, scheme, n)).(1:N)

    G = k -> sum(zip(weights(scheme), shells(scheme))) do (w, shell)
        sum(shell) do b
            A = M[k, k+b]
            q = [angle(A[n, n]) + cartesian(b)' * c[n] for n = 1:N]
            R = hcat([A[:, n] * A[n, n]' for n = 1:N]...)
            T = hcat([(A[:, n] / A[n, n]) * q[n] for n = 1:N]...)
            4w * ((R - R') / 2 - (T + T') / 2im) / length(brillouin_zone)
        end
    end

    return map(G, brillouin_zone)
end

"""
Gauge gradient for the ILA. 

"""
function gauge_gradient(U::Gauge, scheme::W90FiniteDifference, brillouin_zone::B, ::Type{TruncatedConvolution}) where {B<:BrillouinZone}
    M = gauge_transform(neighbor_basis_integral(scheme), U)
    N = n_band(M)

    ρ̃ = Dict(b => (n -> 1 / N * sum(k -> M[k, k+b][n, n], brillouin_zone)).(1:N) for b in vcat(shells(scheme)...))

    G = k -> sum(zip(weights(scheme), shells(scheme))) do (w, shell)
        sum(shell) do b
            R = hcat([M[k, k+b][:, n] * ρ̃[b][n]' / abs(ρ̃[b][n]) for n = 1:N]...)
            2w * (R - R') / length(brillouin_zone)
        end
    end

    return map(G, brillouin_zone)
end


function all_spread(U, scheme, ::Type{Branch}) where {Branch}
    M = gauge_transform(neighbor_basis_integral(scheme), U)
    N = n_band(M)
    return map(1:N) do n
        spread(M, scheme, n, Branch)
    end
end

function all_center(U, scheme, ::Type{Branch}) where {Branch}
    M = gauge_transform(neighbor_basis_integral(scheme), U)
    N = n_band(M)
    return map(1:N) do n
        center(M, scheme, n, Branch)
    end
end

function line_search(U, f, ∇f, ∇f², α; α_0 = 2)
    brillouin_zone = grid(U)
    Ω = f(U)
    ∇Ω = ∇f(U)
    ∇Ω² = ∇f²(∇Ω)

    # println("Ω: $(Ω)"); println("∇Ω²: $(∇Ω²)"); println("α: $(α)");
    # println()

    function try_step(α)
        Q = Ω - 0.5α * ∇Ω²
        new_elements = elements(map(k -> U[k] * cis(-Hermitian(1im * α * ∇Ω[k])), brillouin_zone))
        V = Gauge{typeof(brillouin_zone)}(brillouin_zone, new_elements)
        Ω_V = f(V)
        return V, Ω_V, Q
    end

    while true
        V, Ω_V, Q = try_step(α)

        discontinuous = ∇Ω² > 1e-7 && α < 1e-3
        discontinuous && return let α = α_0, (V, Ω_V, _) = try_step(α)
            println("Discontinuous!")
            println()
            V, Ω_V, α, ∇Ω²
        end
        # println("Q: $(Q)"); println("Ω: $(Ω)"); println("Ω_V: $(Ω_V)"); println("α: $(α)"); flush(stdout)
        # println()
        Ω_V > Q || return V, Ω_V, α, ∇Ω²
        α /= 2
    end
end

struct ILAOptimizer
    scheme::W90FiniteDifference
    meta::Dict
    ILAOptimizer(scheme::W90FiniteDifference) = new(scheme, Dict())
end
scheme(optimizer::ILAOptimizer) = optimizer.scheme

function (optimizer::ILAOptimizer)(U_0::Gauge, ::Type{Branch}; n_iteration = Inf, α_0 = 2) where {Branch}
    N = optimizer |> scheme |> neighbor_basis_integral |> n_band
    U = U_0
    brillouin_zone = grid(U)
    f = U -> sum(all_spread(U, scheme(optimizer), Branch))
    ∇f = U -> Gauge{typeof(brillouin_zone)}(brillouin_zone, elements(gauge_gradient(U, scheme(optimizer), brillouin_zone, Branch)))
    ∇f² = ∇f -> sum(brillouin_zone) do k
        v = reshape(∇f[k], N^2)
        abs(v' * v)
    end
    α = α_0

    current_iteration = -1
    while n_iteration > current_iteration
        current_iteration += 1
        U, Ω, α, ∇Ω² = line_search(U, f, ∇f, ∇f², 2α, α_0 = α_0)
        ∇Ω² < 1e-7 && break
        log(optimizer, U, Branch, current_iteration, Ω, ∇Ω², α)
    end
    println(current_iteration)
    return U
end

function log(optimizer, U, Branch, current_iteration, Ω, ∇Ω², α)
    append!(optimizer.meta[:ila_spread], [all_spread(U, scheme(optimizer), Branch)])
    mod(current_iteration, 10) == 0 || return

    println("Ω: $(Ω)")
    println("∇Ω²: $(∇Ω²)")
    println("α: $(α)")
    u = set_gauge(optimizer.meta[:u], U)
    wannier_orbitals = u(:, optimizer.meta[:phase])
    densities = abs2.(wannier_orbitals)
    # σ_total = 0
    # println("Convolutional: $(σ_total)")
    # @printf "ILA        : %.3f  %.3f  %.3f  %.3f\n" ...
    # @printf "Convolution: "
    center_difference = zeros(N)
    convolutional_spread = zeros(N)
    ila_centers = all_center(U, scheme(optimizer), Branch)
    for i in 1:length(densities)
        ρ = densities[i]
        c, σ = center_spread(fft(ρ, false), optimizer.meta[:r̃2])
        center_difference[i] = norm(c - ila_centers[i])
        convolutional_spread[i] = σ
        # @printf "%.3f  " σ
        # σ_total += σ
    end

    append!(optimizer.meta[:center_difference], [center_difference])
    append!(optimizer.meta[:convolutional_spread], [convolutional_spread])
    println()
end


# if n_iteration == 0
#     for α = -1e-1:5e-3:1e-1
#         _, _, omega, _ = try_step(α)
#         @printf "%.3f  " omega
#     end
#     println()
# end
