export Wannier,
    init_wannier,
    Gauge,
    gauge,
    gauge!,
    set_gauge,
    NeighborIntegral,
    find_neighbors,
    integrals,
    orbital_grid,
    gauge_transform, 
    phase_factors

"""
The gauge ``U^{k}`` is set of matrices, where each matrix
corresponds to each k-point.
"""
struct Gauge <: OnGrid{BrillouinZone}
    grid::BrillouinZone
    elements::Array{AbstractMatrix{ComplexFxx}, <:Any}
end

ket(g::Gauge) = true

Gauge(grid::T) where T <: BrillouinZone = 
    Gauge(grid, Array{Matrix{ComplexFxx}, n_dims(T)}(undef, size(grid)))

function Gauge(grid::T, n::Integer) where T <: BrillouinZone 
    U = Gauge(grid)
    for v in grid
        U[v] = diagm(ones(ComplexFxx, n))
    end
    return U
end

Base.getindex(g::Gauge, k::KPoint) = invoke(getindex, Tuple{OnGrid, KPoint}, g, reset_overflow(k))
Base.setindex!(g::Gauge, value, k::KPoint) = invoke(setindex!, Tuple{OnGrid, Any, KPoint}, g, value, reset_overflow(k))

"""
Wannier is the collection of ``u_{nk}`` orbitals together with a gauge. 

Indexing a wannier with a k-point gives a set of basis functions for the bands at that k-point.
"""
struct Wannier{T <: OnGrid} <: OnGrid{BrillouinZone}
    grid::BrillouinZone
    elements::Array{<:AbstractVector{T}, <:Any}
    gauge::Gauge
end

gauge(wannier::Wannier) = wannier.gauge
gauge!(wannier::Wannier, new_gauge::Gauge) = wannier.gauge = new_gauge
set_gauge(wannier::Wannier, new_gauge::Gauge) = @set wannier.gauge = new_gauge
orbital_grid(wannier::Wannier) = grid(elements(wannier)[1,1,1][1])

function init_wannier(grid::BrillouinZone)
    elements = Array{Vector{UnkBasisOrbital{ReciprocalLattice3D}},3}(undef, size(grid))
    return Wannier(grid, elements, Gauge(grid))
end

function fft(wannier::Wannier{UnkBasisOrbital{T}}) where T <: HomeCell
    g = grid(wannier)
    elements = Array{Vector{UnkBasisOrbital{dual_grid(T)}},n_dims(T)}(undef, size(g))
    transformed = Wannier(g, elements, gauge(wannier))
    Threads.@threads for k in collect(g)
        transformed[k] = fft.(wannier[k])
    end
    return transformed
end

function _fft!(wannier::Wannier{UnkBasisOrbital{T}}) where T <: HomeCell
    Wannier(grid(wannier), 
    elements(map(k->_fft!.(wannier[k]), grid(wannier))),
    gauge(wannier))
end

function ifft(wannier::Wannier{UnkBasisOrbital{T}}) where T <: ReciprocalLattice
    g = grid(wannier)
    elements = Array{Vector{UnkBasisOrbital{dual_grid(T)}},n_dims(T)}(undef, size(g))
    transformed = Wannier(g, elements, gauge(wannier))
    Threads.@threads for k in collect(g)
        transformed[k] = ifft.(wannier[k])
    end
    return transformed
end

function _ifft!(wannier::Wannier{UnkBasisOrbital{T}}) where T <: ReciprocalLattice
    Wannier(grid(wannier), 
    elements(map(k->_ifft!.(wannier[k]), grid(wannier))),
    gauge(wannier))
end

"""
    outer_orbital(unk, G)

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
    return UnkOrbital(g[:, n], wannier[k], true, true)
end

"""
    phase_factors(wannier)

Compute the phase factors eⁱᵏˣ for each k on a limited region of the supercell.
The region is constructing by expanding the homecell by `factors` in each direction.
"""
function phase_factors(wannier::Wannier)
    brillouin_zone = grid(wannier)
    factors = [s for s in size(brillouin_zone)]
    g = orbital_grid(wannier)
    homecell = isa(g, HomeCell) ? g : transform_grid(g)
    supercell = expand(homecell, factors)
    return map(brillouin_zone) do k 
        map(supercell) do r
            exp(1im * k' * r)
        end
    end
end

function (wannier::Wannier)(n::Integer, phase=nothing)
    phase === nothing && (phase = phase_factors(wannier))
    brillouin_zone = grid(wannier)
    N = length(brillouin_zone)

    return (1/N) * sum(brillouin_zone) do k
        phase[k] * expand(wannier[n, k] |> UnkBasisOrbital, [s for s in size(brillouin_zone)])
    end
end

function (wannier::Wannier)(I::Range, phase=nothing)
    phase === nothing && (phase = phase_factors(wannier))
    return [wannier(n, phase) for n in I]
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

"""
The neighbor integrals indexed by two kpoints.
For each pair of kpoint, the integrals are stored 
as a matrix. This is the same as the ``M_{mn}^{k, b}``
matrix in MLWF. The matrix elements are accessed by `M[k, k+b][m, n]`
"""
struct NeighborIntegral
    integrals::Dict{Pair{KPoint,KPoint},Matrix{ComplexFxx}}
end

NeighborIntegral() = NeighborIntegral(Dict())
integrals(n::NeighborIntegral) = n.integrals

function Base.hash(p::Pair{<:KPoint,<:KPoint}) 
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

function Base.setindex!(neighbor_integral::NeighborIntegral, value::AbstractMatrix{ComplexFxx}, g::Vararg{<:KPoint})
    g_1, g_2 = g
    i = integrals(neighbor_integral)
    if haskey(i, g_2 => g_1)
        i[g_2 => g_1] = adjoint(value)
    else
        i[g_1 => g_2] = value
    end
end

"""
    gauge_transform(neighbor_integral, gauge)

Perform a gauge transform on the neighbor integrals.

``U^{k \\dagger} M^{k, k+b} U^{k+b}``

"""
function gauge_transform(neighbor_integral::NeighborIntegral, gauge::Gauge)
    t = NeighborIntegral()
    for ((k_1, k_2), integral) in integrals(neighbor_integral)
        t[k_1, k_2] = adjoint(gauge[k_1]) * integral * gauge[k_2] 
    end
    return t
end
