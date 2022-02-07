export kpoint,
    kpoint!,
    i_kpoint,
    i_kpoint!,
    cache,
    cache!,
    index_band,
    compare,
    AbstractUnkOrbital,
    KPoint,
    UnkBasisOrbital,
    UnkOrbital,
    orthonormal,
    vectorize,
    dagger,
    dagger!

"""
An AbstractUnkOrbital represents the``u_{n, \\mathbf{k}}`` orbitals.

Such a function is also a basis vector and a linear combination.
This is where the single subtyping system of Julia gives me headaches.
"""
abstract type AbstractUnkOrbital{T} <: OnGrid{T} end
const KPoint = AbstractGridVector{<:BrillouinZone}

Base.adjoint(orbital::AbstractUnkOrbital) = dagger(orbital)
i_kpoint(orbital::AbstractUnkOrbital) = haskey(orbital.meta, :i_kpoint) ? orbital.meta[:i_kpoint] : nothing
i_kpoint!(orbital::AbstractUnkOrbital, i) = orbital.meta[:i_kpoint] = i
cache(orbital::AbstractUnkOrbital) = haskey(orbital.meta, :cache) ? orbital.meta[:cache] : nothing
cache!(orbital::AbstractUnkOrbital, c) = orbital.meta[:cache] = c

mutable struct UnkBasisOrbital{T} <: AbstractUnkOrbital{T}
    grid::T
    elements::AbstractArray
    kpoint::AbstractGridVector{<:BrillouinZone}
    index_band::Integer
    ket::Bool
    meta::Dict{Symbol,Any}
end

"""
    UnkBasisOrbital(grid, elements, kpoint, index_band)

This is the concrete type that numerically stores the orbital as a function on a
grid (either the homecell or the reciprocal lattice). It is also associated with
a k-point and an index of a band.

More likely than not, one get these from reading the output of quantum espresso.

```jldoctest orbital
julia> lattice = make_grid(ReciprocalLattice3D, CARTESIAN_BASIS, size_to_domain((4, 4, 4)));

julia> ϕ = map(g->(g==lattice[1, 0, 0]) |> Float64, lattice);

julia> brillouin_zone = make_grid(BrillouinZone3D, CARTESIAN_BASIS, size_to_domain((3, 3, 3)));

julia> ϕ = UnkBasisOrbital(lattice, elements(ϕ), brillouin_zone[1, 0, -1], 2)
ket
grid:
    type: HomeCell3D
    domain: ((-2, 1), (-2, 1), (-2, 1))
    basis:
        ket: 1.000, 0.000, 0.000
        ket: 0.000, 1.000, 0.000
        ket: 0.000, 0.000, 1.000
kpoint:
    GridVector{BrillouinZone3D}:
        coefficients: [1, 0, -1]
    
band:
    2
```
"""
function UnkBasisOrbital(
    grid::T,
    elements::AbstractArray,
    kpoint::AbstractGridVector{<:BrillouinZone},
    index_band::Integer,
) where T 
    orbital = UnkBasisOrbital{T}(grid, elements, kpoint, index_band, true, Dict()) 
    return orbital
end 

"""
    kpoint(orbital)

Get the `kpoint` associated with the `orbital`.

```jldoctest orbital 
julia> kpoint(ϕ)
GridVector{BrillouinZone3D}:
    coefficients: [1, 0, -1]
```
"""
kpoint(orbital::UnkBasisOrbital) = orbital.kpoint

"""
    kpoint!(orbital, new_kpoint)

Set the kpoint associated with `orbital` to `new_kpoint`.

```jldoctest orbital 
julia> kpoint!(ϕ, brillouin_zone[0, 0, -1])
ket
grid:
    type: HomeCell3D
    domain: ((-2, 1), (-2, 1), (-2, 1))
    basis:
        ket: 1.000, 0.000, 0.000
        ket: 0.000, 1.000, 0.000
        ket: 0.000, 0.000, 1.000
kpoint:
    GridVector{BrillouinZone3D}:
        coefficients: [0, 0, -1]
    
band:
    2
```
"""
function kpoint!(orbital::UnkBasisOrbital, new_kpoint::KPoint) 
     orbital.kpoint = new_kpoint
     return orbital
end

"""
    index_band(orbital)

Get the index of the band of the `orbital`.

```jldoctest orbital
julia> index_band(ϕ)
2
```
"""
index_band(orbital::UnkBasisOrbital) = orbital.index_band



## Override the functions of Basis.

"""
    dagger(orbital)

Complex conjutate an orbital. 

```jldoctest orbital
julia> dagger(ϕ) |> ket
false
```
"""
function dagger(orbital::UnkBasisOrbital)
    orbital = @set orbital.ket = !orbital.ket
    orbital = @set orbital.elements = conj.(elements(orbital))
    return orbital
end

"""
    dagger!(orbital)

Complex conjutate an orbital in place. 
"""
function dagger!(orbital::UnkBasisOrbital)
    orbital.elements = conj.(elements(orbital))
    orbital.ket = !orbital.ket
end


function resemble(orbital::UnkBasisOrbital{S}, ::Type{T}, new_elements=nothing) where {S <: Grid,  T <: Grid}
    g = grid(orbital)
    S == dual_grid(T) && (g = transform_grid(g))
    if new_elements === nothing 
       new_elements = zeros(eltype(elements(orbital)), size(g))
    end
    UnkBasisOrbital(g, new_elements, kpoint(orbital), index_band(orbital)) 
end


function Base.show(io::IO, orbital::UnkBasisOrbital)
    ket(orbital) ? print(io, "ket\n") : print(io, "bra\n")
    print(io, "grid:\n$(indent(string(grid(orbital))))" * "\n")
    print(io, "kpoint:\n$(indent(repr(kpoint(orbital))))\n")
    print(io, "band:\n$(indent(repr(index_band(orbital))))")
end

"""
Given a few basis orbitals, we can construct a linear combination of them for a
basic set of symbolic manipulation.  A linear combination is named `UnkOrbital`
due to my stupidity.

This also should be a subtype of OnGrid and Basis, but we are crippled by the
single subtyping system.
"""
mutable struct UnkOrbital <: LinearCombination
    _coefficients::Vector{Number}
    basis::Any
    ket::Bool
    orthonormal::Bool
end

"""
    kpoint(orbital)

Get the k-point.
"""
kpoint(orbital::UnkOrbital) =
    isempty(orbital.basis) ? nothing : kpoint(orbital.basis[1])

"""
    kpoint(orbital)

Set the k-point recursively through the basis vectors.
"""
function kpoint!(orbital::UnkOrbital, new_kpoint::KPoint)
    (b -> kpoint!(b, new_kpoint)).(orbital.basis)
end

"""
    index_band(orbital)

Get the index of the band of the orbital.
"""
index_band(orbital::UnkOrbital) =
    isempty(orbital.basis) ? nothing : index_band(orbital.basis[1])

"""
    orthonormal(orbital)

Whether the basis set is orthonormal.
"""
orthonormal(orbital::UnkOrbital) = orbital.orthonormal

"""
    UnkOrbital(orbital, orthonormal = true)

Construct a linear combination of orbitals from `orbital` as a basis vector.
Set `orthonormal`  to `true` for an orthogonal basis set.

```jldoctest orbital
julia> ϕ₁ = UnkOrbital(map(g->Float64(g==lattice[1, 0, 0]), lattice), true)
ket
coefficients:
    Number[1.0]
n_basis:
    1

julia> ϕ₂ = UnkOrbital(map(g->Float64(g==lattice[0, 1, 0]), lattice), true);

julia> ϕ₃ = UnkOrbital(map(g->Float64(g==lattice[0, 0, 1]), lattice), true); 
```
"""
UnkOrbital(orbital, orthonormal = true) = UnkOrbital(
    [1.0],
    [ket(orbital) ? orbital : dagger(orbital)],
    ket(orbital),
    orthonormal,
)


"""
    dagger(orbital)

Take the complex conjugate recursively through the basis vectors.

```jldoctest orbital
julia> dagger(ϕ₁) |> ket
false
```
"""
function dagger(orbital::UnkOrbital)
    UnkOrbital(
        conj.(coefficients(orbital)),
        _basis(orbital),
        !ket(orbital),
        orthonormal(orbital),
    )
end

"""
    dagger!(orbital)

Take the complex conjugate recursively through the basis vectors in place.
"""
function dagger!(orbital::UnkOrbital)
    orbital._coefficients = conj.(coefficients(orbital))
    orbital.ket = !orbital.ket
end

"""
    braket(o_1, o_2)

Recursively compute the inner product. `o_1` must be a bra and `o_2` must be a
ket (this may be relaxed later). Can also write `o_1 * o_2`

Example:

```jldoctest orbital
julia> ϕ_sum' * ϕ₁
1.0
julia> ϕ_sum' * ϕ₃
0.0
```
"""
function braket(o_1::UnkOrbital, o_2::UnkOrbital)
    !ket(o_1) && ket(o_2) || error("braket requires a bra and a ket.")
    o_1.basis === o_2.basis && orthonormal(o_1) && return transpose(coefficients(o_1)) * coefficients(o_2)

    return invoke(braket, Tuple{LinearCombination,LinearCombination}, o_1, o_2)
end

braket(o_1::UnkBasisOrbital, o_2::UnkOrbital) = braket(UnkOrbital(o_1), o_2)
braket(o_1::UnkOrbital, o_2::UnkBasisOrbital) = braket(o_1, UnkOrbital(o_2))



"""
    add(o_1, o_2, [mutual_orthonormal = false])

Can also write `o_1 + o_2`, where `mutual_orthonormal = false` is assumed.

Add two linear combination of orbitals. The resulting linear combination 
has a orthogonal basis set if both linear combinations have orthogonal 
basis sets and they are mutually orthogonal as specified in `mutual_orthonormal`.

```jldoctest orbital
julia> ϕ_sum = ϕ₁ + ϕ₂
ket
coefficients:
    Number[1.0, 1.0]
n_basis:
    2
julia> orthonormal(ϕ_sum)
false
julia> ϕ_sum = add(ϕ₁, ϕ₂, true)
ket
coefficients:
    Number[1.0, 1.0]
n_basis:
    2

julia> orthonormal(ϕ_sum)
true
```
"""
function add(o_1::UnkOrbital, o_2::UnkOrbital, mutual_orthonormal = false)
    ket(o_1) == ket(o_2) || error("adding a bra to a ket.")
    o_1.basis === o_2.basis &&
        return UnkOrbital(coefficients(o_1) + coefficients(o_2), _basis(o_1), ket(o_1), orthonormal(o_1))

    coefficients_dict = OrderedDict()
    for (k, v) in zip([_basis(o_1)..., _basis(o_2)...], [coefficients(o_1)..., coefficients(o_2)...])
        coefficients_dict[k] = haskey(coefficients_dict, k) ? coefficients_dict[k] + v : v
    end

    # Basis orthogonality is not checked.

    UnkOrbital(
        collect(values(coefficients_dict)),
        collect(keys(coefficients_dict)),
        ket(o_1),
        orthonormal(o_1) && orthonormal(o_2) && mutual_orthonormal,
    )
end

"""
    negate(orbital)

Negate the coefficients of the linear combination. Can also write `-orbital`.
"""
negate(orbital::UnkOrbital) = @set orbital._coefficients = -orbital._coefficients

"""
    mul(s, orbital)

Scale coefficients of the linear combination. Can also write `s * orbital` or `orbital * s`.
"""
mul(s::Number, orbital::UnkOrbital) =
    @set orbital._coefficients = s * orbital._coefficients

"""
    zeros(UnkOrbital, dims...)

Create an array of empty linear combination.  This enables semi-symbolic linear
algebra.

```jldoctest orbital
julia> [ϕ₁, ϕ₂, ϕ₃]' * [0 0 1;
                        0 1 0;
                        1 0 0] * [ϕ₁, ϕ₂, ϕ₃]
1.0
```
"""
function LinearAlgebra.zeros(UnkOrbital, dims...)
    fill(UnkOrbital([], [], true, true), dims...)
end


"""
    zero(UnkOrbital)

Create an empty linear combination. Useful for symbolic linear algebra.
"""
function LinearAlgebra.zero(::UnkOrbital)
    UnkOrbital([], [], true, true)
end

function Base.show(io::IO, orbital::UnkOrbital)
    ket(orbital) ? print(io, "ket\n") : print(io, "bra\n")
    print(io, "coefficients:\n$(indent(repr(coefficients(orbital))))" * "\n")
    print(io, "n_basis:\n$(indent(string(length(basis(orbital)))))" * "\n")
end

function compare(o_1, o_2)
    c_1 = coefficients(kpoint(o_1))
    c_2 = coefficients(kpoint(o_2))
    c_1 == c_2 && return index_band(o_1) < index_band(o_2)
    c_1[1] == c_2[1] && c_1[2] == c_2[2] && return c_1[3] < c_2[3]
    c_1[1] == c_2[1] && return c_1[2] < c_2[2]
    return c_1[1] < c_2[1]
end

consolidate(orbital::UnkBasisOrbital) = (grid(orbital), elements(orbital))

function consolidate(orbital::UnkOrbital)
    elements = nothing
    grid = nothing
    for (c, b) in zip(coefficients(orbital), basis(orbital))
        c == 0 && continue
        new_grid, new_elements = consolidate(b)    
        grid === nothing && (grid = new_grid)
        grid == new_grid || error("Grids Mismatching.")
        elements === nothing && (elements = zeros(ComplexFxx, size(new_grid)))
        elements += c * new_elements
    end
    return grid, elements
end

function UnkBasisOrbital(orbital::UnkOrbital)
    grid, elements = consolidate(orbital)
    UnkBasisOrbital( grid, elements, kpoint(orbital), index_band(orbital))
end
