export OrbitalSet,
    init_orbital_set,
    Gauge,
    random_gauge,
    gauge,
    gauge!,
    set_gauge,
    commit_gauge,
    commit_gauge!,
    orbital_grid,
    gauge_transform,
    phase_factors,
    n_band,
    wannier_orbitals,
    supercell,
    superlattice,
    reciprocal_densities

"""
The gauge ``U^{k}`` is a set of matrices, where each matrix
corresponds to each k-point. The gauge is an instance of `OnGrid`.

By convention (from MLWF), each column of the gauge corresponds to the
coefficients of a transformed ``\\tilde{u}_{n, k}`` orbital. 

row: bands
col: wanniers

`` | \\tilde{u}_{n, k} \\rangle = \\sum_{m} | u_{m, k} \\rangle U^{k}_{m, n} ``
"""
struct Gauge{T<:BrillouinZone} <: OnGrid{T}
    grid::T
    elements::Array{AbstractMatrix{ComplexFxx},<:Any}
end

"""
`OrbitalSet` is the set of ``u_{n, \\mathbf{k}}`` orbitals (`UnkBasisOrbital`) together with a gauge. 
Organizing the ``u_{n, \\mathbf{k}}`` orbitals is very easy with this data structure.
"""
struct OrbitalSet{T<:OnGrid} <: OnGrid{BrillouinZone}
    grid::BrillouinZone
    elements::Array{<:AbstractVector{T},<:Any}
    gauge::Gauge
end

ket(g::Gauge) = true

Gauge(grid::T) where {T<:BrillouinZone} =
    Gauge{T}(grid, Array{Matrix{ComplexFxx},n_dims(T)}(undef, size(grid)))

"""
    Gauge(on_grid)

Construct a gauge from an object on a grid.
"""
Gauge(on_grid::OnGrid{T}) where {T <: BrillouinZone} =
    Gauge{T}(grid(on_grid), elements(on_grid))

"""
    Gauge(brillouin_zone, n)

Create an identity gauge with ``n \\times n`` matrices as the gauge.

Example: 

```jldoctest orbital_set
julia> brillouin_zone = make_grid(BrillouinZone3D, CARTESIAN_BASIS, size_to_domain((4, 4, 4)));

julia> U = Gauge(brillouin_zone, 4)

julia> U[brillouin_zone[0, 0, 1]]
4×4 Matrix{ComplexF64}:
 1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  1.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  1.0+0.0im
```
"""
function Gauge(grid::T, n::Integer) where {T<:BrillouinZone}
    U = Gauge(grid)
    reset_gauge!(U, n)
    return U
end

"""
    random_gauge(brillouin_zone::BrillouinZone, n::Int, ϵ=2π)

Create a random gauge.
"""
function random_gauge(brillouin_zone::BrillouinZone, n::Int, ϵ=2π)
    gauge = Gauge(brillouin_zone, n)
    for k in brillouin_zone
        m = rand(ComplexF64, n, n) * ϵ
        U, _, Vt = svd(m) 
        gauge[k] = U * Vt
    end
    return gauge
end

"""
    reset_gauge!(U::Gauge, n::Integer)

Reset the gauge to ``n \\times n`` identity matrices.
"""
function reset_gauge!(U::Gauge, n::Integer)
    for v in grid(U)
        U[v] = diagm(ones(ComplexFxx, n))
    end
end

function resemble(on_grid::Gauge{S}, ::Type{T}, new_elements = nothing) where {S<:Grid,T<:Grid}
    g = grid(on_grid)
    if new_elements === nothing
        new_elements = zeros(eltype(elements(on_grid)), size(g))
    end
    Gauge{S}(g, new_elements)
end

Base.getindex(g::Gauge, k::KPoint) = invoke(getindex, Tuple{OnGrid,KPoint}, g, reset_overflow(k))
Base.setindex!(g::Gauge, value, k::KPoint) = invoke(setindex!, Tuple{OnGrid,Any,KPoint}, g, value, reset_overflow(k))

"""
    commit_gauge(u)

This will numerically perform the gauge transform on the orbitals in `u`. After the 
orbitals are update, the gauge will be reset to the identity gauge.
"""
commit_gauge(u::OrbitalSet) = commit_gauge!(deepcopy(u))

function commit_gauge!(u::OrbitalSet)
    brillouin_zone = grid(u)
    g = orbital_grid(u)
    tmp = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    function transform!(k)
        Ψ = hcat(vectorize.(u[k])...)
        Ψ = Ψ * gauge(u)[k]
        for (n, c) in enumerate(eachcol(Ψ))
            elements!(u[k][n], reshape(c, size(g)))
        end
        # return [UnkBasisOrbital(reciprocal_lattice, reshape(c, size(reciprocal_lattice)...), k, n) for (n, c) in enumerate(eachcol(Ψ))]
    end

    for k in collect(brillouin_zone)
        transform!(k)
    end
    reset_gauge!(gauge(u), n_band(u))
    BLAS.set_num_threads(tmp)
    return u
end


"""
    gauge(orbital_set)

Get the gauge of an orbital set.
"""
gauge(orbital_set::OrbitalSet) = orbital_set.gauge
gauge!(orbital_set::OrbitalSet, new_gauge::Gauge) = orbital_set.gauge = new_gauge

"""
    set_gauge(orbital_set, new_gauge)
    
Create a new orbital set with the gauge replaced by `new_gauge`.
"""
set_gauge(orbital_set::OrbitalSet, new_gauge::Gauge) = @set orbital_set.gauge = new_gauge


function init_orbital_set(grid::BrillouinZone)
    elements = Array{Vector{UnkBasisOrbital{ReciprocalLattice3D}},3}(undef, size(grid))
    return OrbitalSet(grid, elements, Gauge(grid))
end

"""
    init_orbital_set(brillouin_zone, ::Type{T})

Create an empty orbital set where the k-points are in the `brillouin_zone` and 
the orbitals are of type `T`. The default type is `UnkBasisOrbital{ReciprocalLattice3D}` 
is you omit the type.

Example:

```jldoctest orbital_set
julia> u = init_orbital_set(brillouin_zone);

julia> typeof(u)
OrbitalSet{UnkBasisOrbital{ReciprocalLattice3D}}
```

If you are reading orbitals from Quantum Espresso, you should do the following instead.
```jldoctest orbital_set
path_to_si = "../test/test_data/test_5"
wave_functions_list = wave_functions_from_directory(joinpath(path_to_si, "si.save"))
k_map, brillouin_zone = i_kpoint_map(wave_functions_list)
ũ = orbital_set_from_save(wave_functions_list)
ũ[brillouin_zone[0, 0, 0]][1]

# output

ket
grid:
    type: ReciprocalLattice3D
    domain: ((-12, 11), (-12, 11), (-12, 11))
    basis:
        ket: -0.612, -0.612, 0.612
        ket: 0.612, 0.612, 0.612
        ket: -0.612, 0.612, -0.612
kpoint:
    GridVector{BrillouinZone3D}:
        coefficients: [0, 0, 0]
    
band:
    1

```
"""
function init_orbital_set(grid::S, ::Type{T}) where {S <: Grid, T<:OnGrid}
    elements = Array{Vector{T}, n_dims(S)}(undef, size(grid))
    return OrbitalSet(grid, elements, Gauge(grid))
end

"""
    orbital_grid(orbital_set)

Get the common grid of the orbitals contained by the `orbital_set`.

```jldoctest orbital_set
julia> lattice = orbital_grid(ũ)
type: ReciprocalLattice3D
domain: ((-12, 11), (-12, 11), (-12, 11))
basis:
    ket: -0.612, -0.612, 0.612
    ket: 0.612, 0.612, 0.612
    ket: -0.612, 0.612, -0.612
```
"""
orbital_grid(orbital_set::OrbitalSet) = grid(elements(orbital_set)[1, 1, 1][1])

"""
    superlattice(ũ)

Get the superlattice associated with `ũ`, which is the brillouin_zone copied 
over the entire crystal lattice.
"""
superlattice(orbital_set::OrbitalSet{UnkBasisOrbital{T}}) where T <:ReciprocalLattice=
expand(grid(orbital_set), [size(orbital_grid(orbital_set))...])

superlattice(orbital_set::OrbitalSet{UnkBasisOrbital{T}}) where T <: HomeCell =
expand(transform_grid(grid(orbital_set)), [size(orbital_grid(orbital_set))...])

"""
    supercell(u)

Get the supercell associated with `u`, which is the homecell copied 
over the entire crystal lattice.
"""
supercell(orbital_set::OrbitalSet{UnkBasisOrbital{T}}) where T <: HomeCell =
expand(orbital_grid(orbital_set), [size(grid(orbital_set))...])

supercell(orbital_set::OrbitalSet{UnkBasisOrbital{T}}) where T <: ReciprocalLattice =
expand(transform_grid(orbital_grid(orbital_set)), [size(grid(orbital_set))...])

"""
    n_band(orbital_set)

Get the number of bands in the set.

```jldoctest orbital_set
julia> n_band(ũ)
4
```
"""
n_band(orbital_set::OrbitalSet) = length(elements(orbital_set)[1, 1, 1])

"""
    fft(u)

Perform an fft on every orbital in the set.
"""
function fft(orbital_set::OrbitalSet{UnkBasisOrbital{T}}) where {T<:HomeCell}
    g = grid(orbital_set)
    elements = Array{Vector{UnkBasisOrbital{dual_grid(T)}},n_dims(T)}(undef, size(g))
    transformed = OrbitalSet(g, elements, gauge(orbital_set))
    Threads.@threads for k in collect(g)
        transformed[k] = fft.(orbital_set[k])
    end
    return transformed
end

function _fft!(wannier::OrbitalSet{UnkBasisOrbital{T}}) where {T<:HomeCell}
    OrbitalSet(grid(wannier),
        elements(map(k -> _fft!.(wannier[k]), grid(wannier))),
        gauge(wannier))
end


"""
    fft(ũ)

Perform an inverse fft on every orbital in the set.

Example: 

```jldoctest orbital_set
julia> u = ifft(ũ)
julia> orbital_grid(u)
```
"""
function ifft(wannier::OrbitalSet{UnkBasisOrbital{T}}) where {T<:ReciprocalLattice}
    g = grid(wannier)
    elements = Array{Vector{UnkBasisOrbital{dual_grid(T)}},n_dims(T)}(undef, size(g))
    transformed = OrbitalSet(g, elements, gauge(wannier))
    Threads.@threads for k in collect(g)
        transformed[k] = ifft.(wannier[k])
    end
    return transformed
end

function _ifft!(wannier::OrbitalSet{UnkBasisOrbital{T}}) where {T<:ReciprocalLattice}
    OrbitalSet(grid(wannier),
        elements(map(k -> _ifft!.(wannier[k]), grid(wannier))),
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

"""
    phase_factors(wannier)

Compute the phase factors eⁱᵏˣ for each k on a limited region of the supercell.
The region is constructing by expanding the homecell by `factors` in each direction.
"""
function phase_factors(wannier::OrbitalSet)
    brillouin_zone = grid(wannier)
    factors = [size(brillouin_zone)...]
    g = orbital_grid(wannier)
    homecell = isa(g, HomeCell) ? g : transform_grid(g)
    supercell = expand(homecell, factors)
    k_coordinates = hcat(cartesian.(brillouin_zone(1:length(brillouin_zone)))...)
    r_coordinates = hcat(cartesian.(supercell(1:length(supercell)))...)
    @time phase = exp.(r_coordinates' * (1im * k_coordinates))
    @time SimpleFunctionOnGrid(brillouin_zone, reshape((n ->
                SimpleFunctionOnGrid(supercell, reshape(phase[:, n], size(supercell)),
                    true)).(1:length(brillouin_zone)), size(brillouin_zone)), true)

    # return map(brillouin_zone) do k 
    #     map(supercell) do r
    #         exp(1im * k' * r)
    #     end
    # end
end

"""
    getindex(u, k)

Indexing a wannier object with a kpoint gives the set of basis orbitals at that
kpoint. The gauge is irrelevant for this type of indexing. One can also write `u[k]`.

Example: 

```jldoctest orbital_set
julia> typeof(ũ[brillouin_zone[0, 0, 0]])
Vector{UnkBasisOrbital{ReciprocalLattice3D}} (alias for Array{UnkBasisOrbital{ReciprocalLattice3D}, 1})
```

The kpoint can be out of the first Brillouin zone. In that case, the basis
orbital corresponding to that kpoint will be computed by phase shifting its
image in the first Brillouin zone.

```jldoctest orbital_set
julia> k_1, k_2 = brillouin_zone[-2, 0, 0], brillouin_zone[2, 0, 0];
julia> has_overflow(k_2)
true
julia> ũ[k_2][1][lattice[1:4, 0, 0]]
4-element Vector{ComplexF64}:
   -0.01524761578943151 - 0.0422291275292075im
  0.0023520566748192902 - 0.005010400928812839im
   0.001762807696221052 - 0.0006365001497096597im
 -0.0001863994930269115 - 8.750521988524712e-5im

 julia> ũ[k_1][1][lattice[1:4, 0, 0]]
 4-element Vector{ComplexF64}:
    -0.6169403121751041 - 0.28961297828233756im
   -0.01524761578943151 - 0.0422291275292075im
  0.0023520566748192902 - 0.005010400928812839im
   0.001762807696221052 - 0.0006365001497096597im
```
"""
function Base.getindex(u::OrbitalSet, k::KPoint)
    flagged_overflow = has_overflow(k)
    images_in_brillouin_zone = invoke(Base.getindex, Tuple{OnGrid,KPoint}, u, flagged_overflow ? reset_overflow(k) : k)
    translated_orbitals = flagged_overflow ? (o -> o >> -overflow(k)).(images_in_brillouin_zone) : images_in_brillouin_zone
    return (o -> kpoint!(o, k)).(translated_orbitals)
end

"""
    getindex(u, n, k)

This fetches us the ``u_{n, \\mathbf{k}}`` orbitals with the gauge applied.
Can also write `u[n, k]`. The result will be a linear combination of orbitals.

Example:

```jldoctest orbital_set
amn = AMN(joinpath(path_to_si, "output/pw2wan/si.amn"))
U = Gauge(grid(ũ), amn, k_map)
ũ = set_gauge(ũ, U)
ũ[1, brillouin_zone[0, 0, 0]]

# output

ket
coefficients:    
    Number[-0.4240209935568641 + 0.06876931361986612im, 0.007756005591334594 - 0.28204917169307464im, 0.8377575800925453 + 0.1167214564594588im, 0.06944769261399274 + 0.12482164973224028im]
n_basis:    
    4
```
"""
function Base.getindex(orbital_set::OrbitalSet, n::Integer, k::KPoint)
    g = gauge(orbital_set)[k]
    return UnkOrbital(g[:, n], orbital_set[k], true, true)
end


function (u::OrbitalSet)(n::Integer, phase = nothing)
    phase === nothing && (phase = phase_factors(u))
    brillouin_zone = grid(u)
    N = length(brillouin_zone)

    return (1 / N) * sum(brillouin_zone) do k
        phase[k] * expand(u[n, k] |> UnkBasisOrbital, [size(brillouin_zone)...])
    end
end

# function (wannier::OrbitalSet)(I::Range, phase = nothing)
#     phase === nothing && (phase = phase_factors(wannier))
#     return [wannier(n, phase) for n in I]
# end
# if coefficients(grid_vector) == [11, 17, 45] || coefficients(grid_vector) == [18, 17, 45]
#     println("k: $(k)\ng:$(g)")
#     println(transformed[k][n][g])
# end
# function bloch_orbital_sum(n) 
#     (1/N) * sum(brillouin_zone) do k
#         expanded = expand(transformed[k][n], [size(brillouin_zone)...])
#         phase[k] * expanded
#     end
# end


"""
    wannier_orbitals(ũ)

Construct all the wannier orbitals in the reciprocal space by interleaving
the ``u_{n, \\mathbf{k}}`` orbtials. Can also write `ũ(:)`. Gettting them 
one at a time is currently not supported due to pathetic performance.

Keep in mind that this is likely the most costy operation in this package even
though the implementation is linear time.

Example: 

```jldoctest orbital_set
julia> wanniers = ũ(:);
julia> wanniers[1]
ket
grid:
    type: ReciprocalLattice3D
    domain: ((-48, 47), (-48, 47), (-48, 47))
    basis:
        ket: -0.153, -0.153, 0.153
        ket: 0.153, 0.153, 0.153
        ket: -0.153, 0.153, -0.153
kpoint:
    GridVector{BrillouinZone3D}:
        coefficients: [0, 0, 0]
    
band:
    1
```
"""
wannier_orbitals(ũ::OrbitalSet) = ũ(:)

function (ũ::OrbitalSet{UnkBasisOrbital{T}})(::Colon) where {T<:ReciprocalLattice}
    # phase === nothing && (phase = phase_factors(u))
    brillouin_zone = grid(ũ)
    reciprocal_lattice = orbital_grid(ũ)

    # reciprocal_supercell = expand(make_grid(typeof(reciprocal_lattice), basis(brillouin_zone), domain(brillouin_zone)), [size(reciprocal_lattice)...])
    reciprocal_supercell = supercell(ũ) |> transform_grid

    orbital_elements = zeros(ComplexFxx, length(brillouin_zone) * length(reciprocal_lattice), n_band(ũ))
    Threads.@threads for k in collect(brillouin_zone)
    # for k in collect(brillouin_zone)
        U = hcat(vectorize.(ũ[k])...)
        for g in collect(reciprocal_lattice)
            grid_vector = reset_overflow(snap(reciprocal_supercell, cartesian(g) + cartesian(k)))
            orbital_elements[linear_index(grid_vector), :] = U[linear_index(g), :]
        end
    end

    return [UnkBasisOrbital(reciprocal_supercell, reshape(orbital_elements[:, n],
            size(reciprocal_supercell)), brillouin_zone(1), n) |> square_normalize! for n in 1:n_band(ũ)]
end

"""
    reciprocal_densities(ũ)

Gives the reciprocal space densities for the wannier functions.
"""
function reciprocal_densities(ũ::OrbitalSet{UnkBasisOrbital{T}}) where {T<:ReciprocalLattice}
    wanniers = commit_gauge(ũ)(:)
    return (ρ->fft(ρ, false)).(abs2.(ifft.(wanniers)))
end