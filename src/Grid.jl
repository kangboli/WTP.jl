export Grid, domain, domain_matrix, set_domain, basis, dual_grid,
    HomeCell, ReciprocalLattice, BrillouinZone, RealLattice, transform_grid, snap,
    x_min, x_max, y_min, y_max, z_min, z_max,
    mins, maxes, HomeCell3D, ReciprocalLattice3D, BrillouinZone3D, RealLattice3D, n_dims, array, invert_grid, shrink, make_grid

"""
Non-orthogonal 3D periodic grid that comply to the centering convention.

It's made of a set of basis vectors and the domain in units of the basis
vectors.  If you subtype this, the concrete type should have these two fields
for things to work out of the box.

```julia
basis::Tuple{Vector3, Vector3, Vector3}
domain::Tuple{Tuple{Integer, Integer}, Tuple{Integer, Integer}, Tuple{Integer, Integer}}
```
"""
abstract type Grid end

# The Domain

"""
    domain(grid)

The domain of the grid are the extremal miller indices.
The grid includes the end points of the domain.

The result is a tuple of `n_dim` tuples in the form
`((x_min, x_max), (y_min, y_max), (z_min, z_max))`.

"""
domain(grid::Grid) = grid.domain

_convert_domain_to_matrix(domain) = SMatrix{3,2,Int}(vcat([[d[1] d[2]] for d in domain]...))

domain_matrix(grid::Grid) = grid._cache._domain_matrix

"""
    set_domain(grid, new_domain)

Create a copy of `grid` with its domain set to `new_domain`.
"""
function set_domain(grid::Grid, new_domain::Tuple)
    grid = @set grid.domain = new_domain
    size = _compute_size(new_domain)
    domain_matrix = _convert_domain_to_matrix(new_domain)
    grid = @set grid._cache = GridCache(
        basis_matrix(grid),
        domain_matrix,
        size,
        _compute_center(grid._cache._size, new_domain),
        tuple(domain_matrix[:, 1]...),
        tuple(domain_matrix[:, 2]...)
    )
    return grid
end


x_min(grid::Grid) = mins(grid)[1]
x_max(grid::Grid) = maxes(grid)[1]
y_min(grid::Grid) = mins(grid)[2]
y_max(grid::Grid) = maxes(grid)[2]
z_min(grid::Grid) = mins(grid)[3]
z_max(grid::Grid) = maxes(grid)[3]

"""
    expand(grid, factors=[2, 2, 2])

Expand a grid by a factor of `factor`. The returned grid will have the 
save basis vector, but it will be larger.
"""
function expand(grid::Grid, factors::Vector{Int} = [2, 2, 2])
    return set_domain(grid, size_to_domain(factors .* size(grid)))
end

"""
    mins(grid::Grid)

The minimum indices of the grid along each direction.

Example: 

```julia
homecell = make_grid(HomeCell3D, CARTESIAN_BASIS, ((-2, 1), (-2, 1), (-2, 1)))
mins(homecell) # gives [-2, -2, -2]
```
"""
mins(grid::Grid) = grid._cache._mins # domain_matrix(grid)[:, 1]

"""
    maxes(grid)

Similar to `mins(grid)`, but gives the maximum indices instead.
"""
maxes(grid::Grid) = grid._cache._maxes # domain_matrix(grid)[:, 2]

"""
    n_dims(T)

Gives the number of dimensions of the grid.

Example:

```julia
n_dims(HomeCell3D) # 3
```
"""
n_dims(::Type{<:Grid}) = 3

"""
    center(grid)

The center of a grid. This should be a tuple of zeros for a grid that complies 
to the centering convention.  The result will be a tuple of `n_dims` integers.
A tuple instead of a vector is used for performance reasons. 
"""
center(grid::Grid)::Tuple = grid._cache._center
_compute_center(size::Tuple, domain::Tuple)::Tuple = map(_center_1d, size, domain)
_center_1d(s, d)::Int = iseven(s) ? (sum(d) + 1) ÷ 2 : sum(d) ÷ 2

# function center(grid::T) where T <: Grid
#     result = zeros(Int, n_dims(T))
#     s, d = size(grid), domain(grid)
#     for i = 1:n_dims(T)
#         u, l = d[i]
#         result[i] = iseven(s[i]) ? (u + l + 1) ÷ 2 : (u + l) ÷ 2
#     end
#     return result
# end


"""
    basis(g)

Basis vectors for the grid. Depending on the type of grid,
the basis is that of the grid such as real/reciprocal lattice vectors.

A tuple of vectors will be returned, and the vectors returned are always ket
vectors.

If it is a matrix that you want, use `basis_matrix` instead for better
performance.
"""
basis(grid::Grid) = grid.basis

"""
    basis_matrix(g)

The basis vector as a matrix, where each column is a basis vector.
The matrix is cached because this has an impact on performance.
"""
basis_matrix(grid::Grid) = grid._cache._basis_matrix

Base.:(==)(g_1::T, g_2::T) where {T<:Grid} = basis(g_1) == basis(g_2) && domain(g_1) == domain(g_2)

_compute_size(domain) = Tuple(d[2] - d[1] + 1 for d in domain)

"""
    size(g)

The size of the grid. The result is a tuple of `n_dims` integers.
"""
Base.size(g::Grid)::Tuple = g._cache._size

"""
    length(g)

The number of grid points in the grid.
"""
Base.length(g::Grid)::Int = prod(size(g))

const Range = AbstractVector{<:Integer}
const NotColon = Union{Range,Integer}

"""
    g[0, 0, 1]

Indexing a grid with a miller index gives the grid vector corresponding to that
index.
"""
Base.getindex(g::Grid, i::Integer, j::Integer, k::Integer) = make_grid_vector(g, SVector(i, j, k))

Base.getindex(g::Grid, I::Range, J::Range, K::Range) = [g[i, j, k] for i in I, j in J, k in K]

Base.getindex(g::Grid, i::Integer, J::Range, K::Range) = [g[i, j, k] for j in J, k in K]
Base.getindex(g::Grid, I::Range, j::Integer, K::Range) = [g[i, j, k] for i in I, k in K]
Base.getindex(g::Grid, I::Range, J::Range, k::Integer) = [g[i, j, k] for i in I, j in J]

Base.getindex(g::Grid, I::Range, j::Integer, k::Integer) = [g[i, j, k] for i in I]
Base.getindex(g::Grid, i::Integer, J::Range, k::Integer) = [g[i, j, k] for j in J]
Base.getindex(g::Grid, i::Integer, j::Integer, K::Range) = [g[i, j, k] for k in K]

Base.getindex(g::Grid, ::Colon, J::Any, K::Any) = g[x_min(g):x_max(g), J, K]
Base.getindex(g::Grid, I::NotColon, ::Colon, K::Any) = g[I, y_min(g):y_max(g), K]
Base.getindex(g::Grid, I::NotColon, J::NotColon, ::Colon) = g[I, J, z_min(g):z_max(g)]

"""
Indexing a grid with a single integer i gives a grid vector. 
Indexing an OnGrid object with this grid vector gives the ith element
of the object.

>>> elements(o)[i] == o[g[i]]

This provides a mapping between 1D indices (which are to be used for matrix algorithm) 
and grid vectors.
"""
(g::Grid)(i::Integer) = make_grid_vector(g,
    SVector(standard_to_miller(size(g), one_to_three(i, size(g)), (0, 0, 0))))

"""
Indexing a grid with a list of linear indices gives a list of grid vector.
"""
(g::Grid)(linear_indices::Range) = [g(i) for i in linear_indices]
(g::Grid)(::Colon) = [g(i) for i = 1:length(g)]

"""
    array(grid)

Convert the grid to a multidimensional array.


Example:

```julia
homecell = make_grid(HomeCell3D, CARTESIAN_BASIS, ((-2, 1), (-2, 1), (-2, 1)))
array(homecell)[1, 1, 1]
# grid: HomeCell3D
# _coefficients: [-2, -2, -2]
```
"""
array(grid::Grid) = grid[:, :, :]

# Iteration

"""
Iterate over the grid gives a sequence of grid vectors that goes through each
grid point. The order of the iteration is not guaranteed.
"""
function Base.iterate(grid::Grid)
    linear_indices = 1:length(grid)
    first, linear_indices_state = iterate(linear_indices)
    return (grid(first), (linear_indices, linear_indices_state))
end

function Base.iterate(grid::Grid, state)
    linear_indices, linear_indices_state = state
    next_linear_index = iterate(linear_indices, linear_indices_state)
    next_linear_index === nothing && return nothing
    next, linear_indices_state = next_linear_index
    return (grid(next), (linear_indices, linear_indices_state))
end

"""
    dual_grid(T)

The type of the dual grid of T.
"""
dual_grid(::Type{T}) where {T} = T

"""
    transform_grid(g)

FFT of the grid.

The resulting basis vectors should satisfy
aᵢᵀ bᵢ = 2π / sᵢ,
where s is the size of the grid (number of grid points in each direction).
"""
transform_grid(grid::T) where {T<:Grid} =
    let A = basis_matrix(grid) * diagm([size(grid)...]), B = 2 * pi * inv(A)'
        make_grid(dual_grid(T), matrix_to_vector3(B), size_to_domain(size(grid)))
    end

"""
    snap(grid, point)

Snap a coordinate to a grid point.  
"""
snap(grid::Grid, point::AbstractVector) =
    make_grid_vector(grid, (i -> round(Int, i)).((basis_matrix(grid) \ point)))


"""
Two grids are identical if they are the same concrete type; they have the same
domain; and they have the same basis.
"""
Base.:(==)(grid_1::Grid, grid_2::Grid) =
    typeof(grid_1) == typeof(grid_2) &&
    domain(grid_1) == domain(grid_2) &&
    basis(grid_1) == basis(grid_2)


function Base.string(grid::Grid)::String
    type_str = "type: $(typeof(grid))"
    domain_str = "domain: $(domain(grid))"
    grid_basis_str = "basis:\n" * join(["$(indent(repr(b)))" for b in basis(grid)], "\n")
    return join([type_str, domain_str, grid_basis_str], "\n")
end

function html(grid::Grid)
    domain_html = "<ul>" * join(["<li>from $(d[1]) to $(d[2])</li>" for d in domain(grid)], "\n") * "</ul>"
    basis_html = "<ul>" * join(["<li>$(html(b))</li>" for b in basis(grid)], "\n") * "</ul>"
    return "<ul>
    <li>Type: $(typeof(grid))</li>
    <li>Domain: $(domain_html)</li>
    <li>Basis: $(basis_html)</li>
    </ul>"
end

function Base.show(io::IO, ::MIME"text/plain", grid::Grid)
    println(io, string(grid))
end

function Base.show(io::IO, ::MIME"text/html", grid::Grid)
    println(io, html(grid))
end

mutable struct GridCache
    _basis_matrix::SMatrix{3,3,Float64}
    _domain_matrix::SMatrix{3,2,Int}
    _size::NTuple{3,Integer}
    _center::NTuple{3,Integer}
    _mins::NTuple{3,Integer}
    _maxes::NTuple{3,Integer}
end

"""
There are four basic grids (two pairs of dual grids) in Condensed Phase on which
most physical concepts are defined. Generic grid operations should be defined
on the abstract type `Grid`. Operations specific to each space should be defined
on the concrete grid.

This approach is to be contrasted with allocating a Fortran array and program  
functions as array manipulation.
"""

"""
The abstract homecell. It is reciprocal to the reciprocal lattice.
"""
abstract type HomeCell <: Grid end

"""
The home cell in 3D. The ``u_{nk}`` orbitals in the real space is represented as
function on this grid. 
"""
struct HomeCell3D <: HomeCell
    basis::NTuple{3,Vector3}
    domain::NTuple{3,NTuple{2,Integer}}
    _cache::GridCache
end

abstract type ReciprocalLattice <: Grid end

"""
The 3D reciprocal lattice. The ``u_{nk}`` orbitals in the frequency space are
represented as functions on this grid.

This grid is the dual grid of `HomeCell3D`.
Given a `homecell`, one can find its corresponding reciprocal lattice by

```julia
reciprocal_lattice = transform_grid(homecell)
```
"""
struct ReciprocalLattice3D <: ReciprocalLattice
    basis::NTuple{3,Vector3}
    domain::NTuple{3,NTuple{2,Integer}}
    _cache::GridCache
end

abstract type BrillouinZone <: Grid end
"""
The 3D Brillouin zone. The usage is the same of `HomeCell3D`.
"""
struct BrillouinZone3D <: BrillouinZone
    basis::NTuple{3,Vector3}
    domain::NTuple{3,NTuple{2,Integer}}
    _cache::GridCache
end

abstract type RealLattice <: Grid end
"""
The 3D crystal lattice. This is the dual grid of `BrillouinZone3D`.
"""
struct RealLattice3D <: RealLattice
    basis::NTuple{3,Vector3}
    domain::NTuple{3,NTuple{2,Integer}}
    _cache::GridCache
end

function make_grid(::Type{T}, basis, domain) where {T<:Grid}
    size = _compute_size(domain)
    domain_matrix = _convert_domain_to_matrix(domain)
    cache = GridCache(vector3_to_matrix(basis),
        domain_matrix,
        size,
        _compute_center(size, domain),
        tuple(domain_matrix[:, 1]...),
        tuple(domain_matrix[:, 2]...)
    )

    return T(basis, domain, cache)
end

dual_grid(::Type{HomeCell3D}) = ReciprocalLattice3D
dual_grid(::Type{ReciprocalLattice3D}) = HomeCell3D
dual_grid(::Type{BrillouinZone3D}) = RealLattice3D
dual_grid(::Type{RealLattice3D}) = BrillouinZone3D

# inverse_grid(::Type{HomeCell3D}) = RealLattice3D
# inverse_grid(::Type{RealLattice3D}) = HomeCell3D
# inverse_grid(::Type{BrillouinZone3D}) = ReciprocalLattice3D
# inverse_grid(::Type{ReciprocalLattice3D}) = BrillouinZone3D
