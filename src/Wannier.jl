export Wannier,
    init_wannier,
    Gauge,
    gauge,
    NeighborIntegral,
    find_neighbors,
    integrals,
    reciprocal_lattice,
    gauge_transform

struct Gauge <: OnGrid{BrillouinZone}
    grid::BrillouinZone
    elements::Array{AbstractMatrix{ComplexF32}, 3}
end

Gauge(grid::BrillouinZone) = Gauge(grid, Array{Matrix{ComplexF32},3}(undef, size(grid)))

Base.getindex(g::Gauge, k::KPoint) = invoke(getindex, Tuple{OnGrid, KPoint}, g, reset_overflow(k))
Base.setindex!(g::Gauge, value, k::KPoint) = invoke(setindex!, Tuple{OnGrid, Any, KPoint}, g, value, reset_overflow(k))

"""
Wannier represents a Brillouin zone of bands of orbitals.

Indexing a wannier with a k-point gives a set of basis functions 
for the bands at that k-point.
"""
struct Wannier{T <: OnGrid} <: OnGrid{BrillouinZone}
    grid::BrillouinZone
    elements::Array{Vector{T},3}
    gauge::Gauge
end

gauge(wannier::Wannier) = wannier.gauge
reciprocal_lattice(wannier::Wannier) = grid(elements(wannier)[1,1,1][1])

function init_wannier(grid::BrillouinZone)
    elements = Array{Vector{UnkBasisOrbital{ReciprocalLattice}},3}(undef, size(grid))
    return Wannier(grid, elements, Gauge(grid))
end

function fft(wannier::Wannier{UnkBasisOrbital{HomeCell}}) 
    g = grid(wannier)
    elements = Array{Vector{UnkBasisOrbital{dual_grid(HomeCell)}},3}(undef, size(g))
    transformed = Wannier(g, elements, gauge(wannier))
    for k in g
        transformed[k] = fft.(wannier[k])
    end
    return transformed
end

function ifft(wannier::Wannier{UnkBasisOrbital{ReciprocalLattice}}) 
    g = grid(wannier)
    elements = Array{Vector{UnkBasisOrbital{dual_grid(ReciprocalLattice)}},3}(undef, size(g))
    transformed = Wannier(g, elements, gauge(wannier))
    for k in g
        transformed[k] = ifft.(wannier[k])
    end
    return transformed
end


"""
Get an orbital out of the brillouin zone from its copy within the brillouin zone
through a translation in frequency space.

ψ_{nk} = e^{ikx} u_{nk}(x) = ψ_{nk+G} = e^{i(k+G)x} u_{nk+G}(x)
u_{nk+G}(x) = e^{-iGx} u_{nk}(x)

Here, G is a vector in the recirpocal lattice basis.
"""
function outer_orbital(unk::UnkBasisOrbital{<:Grid}, G)::UnkBasisOrbital
    unkg = translate(unk, -grid(unk)[G...])
    kpoint!(unkg, add_overflow(kpoint(unkg), G))
    return unkg
end

function Base.getindex(wannier::Wannier, n::Integer, k::KPoint)
    # gauge = elements(gauge(wannier))[miller_to_standard(k, translation(wannier))...]
    g = gauge(wannier)[k]
    return UnkOrbital(g[n, :], wannier[k], true, true)
end

"""
Indexing a wannier object with a kpoint gives the set of basis 
orbitals at that kpoint.

The kpoint can be out of the first Brillouin zone. In that case, 
the basis orbital corresponding to that kpoint will be computed
by phase shifting its image in the first Brillouin zone.

"""
function Base.getindex(u::Wannier, k::KPoint)
    images_in_brillouin_zone = invoke(Base.getindex, Tuple{OnGrid,KPoint}, u, reset_overflow(k))
    translated_orbitals = has_overflow(k) ?  
        (o -> standardize(translate(o, -grid(o)[overflow(k)...]))).(images_in_brillouin_zone) :
        images_in_brillouin_zone
    return (o -> kpoint!(o, k)).(translated_orbitals)
end

struct NeighborIntegral
    integrals::Dict{Pair{KPoint,KPoint},Matrix{ComplexF32}}
end

NeighborIntegral() = NeighborIntegral(Dict())
integrals(n::NeighborIntegral) = n.integrals

function Base.hash(p::Pair{KPoint,KPoint}) 
    m, n = p
    return hash(m) + hash(n)
end

function Base.getindex(neighbor_integral::NeighborIntegral, k_1::KPoint, k_2::KPoint)
    k_1 == k_2 && return I
    i = integrals(neighbor_integral)
    haskey(i, k_1 => k_2) && return i[k_1 => k_2]
    haskey(i, k_2 => k_1) && return adjoint(i[k_2 => k_1])
    return nothing
end

function Base.setindex!(neighbor_integral::NeighborIntegral, value::AbstractMatrix{ComplexF32}, g::Vararg{KPoint})
    g_1, g_2 = g
    i = integrals(neighbor_integral)
    if haskey(i, g_2 => g_1)
        i[g_2 => g_1] = adjoint(value)
    else
        i[g_1 => g_2] = value
    end
end

function gauge_transform(neighbor_integral::NeighborIntegral, gauge::Gauge)
    t = NeighborIntegral()
    for ((k_1, k_2), integral) in integrals(neighbor_integral)
        t[k_1, k_2] = adjoint(gauge[k_1]) * integral * gauge[k_2] 
    end
    return t
end


    # haskey(neighbor_set, g1) && haskey(neighbor_set[g1], g2) && return neighbor_set[g1][g2]
    # haskey(neighbor_set, g2) &&
    #     haskey(neighbor_set[g2], g1) &&
    #     return adjoint(neighbor_set[g2][g1])

# b_vectors = vcat(neighbor_shells[1:n_shell]...)


    # f(k) = UnkOrbital(wannier[n, k])
    # g(x...) = grid_vector_coefficientsonstructor(grid(k), x)
    # c(b) = basis_transform(coefficients(b), basis(b), CARTESIAN_BASIS)

    # b_1, b_2, b_3 = g(1, 0, 0), g(0, 1, 0), g(0, 0, 1)
    # B = [c(b_1) c(b_2) c(b_3)]'
    # # B = vector3_to_matrix(basis(g))'
    # rhs = 1 / 2 * [f(k + b_1) - f(k - b_1), f(k + b_2) - f(k - b_2), f(k + b_3) - f(k - b_3)]
    # # rhs = [f(k + b_1) - f(k), f(k + b_2) - f(k), f(k + b_3) - f(k)]
    # return inv(B) * rhs


# First order scheme. C₄² = 6 integrals per point.
# | b₁ᵀ | ⋅ ∇f(k) =  | f(k+b₁)-f(k) |
# | b₂ᵀ |            | f(k+b₂)-f(k) |
# | b₃ᵀ |            | f(k+b₃)-f(k) |

# Second order scheme. C₆² = 15 integrals per point.
# | b₁ᵀ | ⋅ ∇f(k) = 1/2 | f(k+b₁)-f(k-b₁) |
# | b₂ᵀ |                | f(k+b₂)-f(k-b₂) |
# | b₃ᵀ |                | f(k+b₃)-f(k-b₃) |

# Wannier90 requires 6-12 integrals and computes only to the first order.


# neighbor_shells = find_shells(wannier, n_shells)
# weights = compute_weights(neighbor_shells)

# weights === nothing && return finite_difference(wannier, n, k, n_shells+1)

