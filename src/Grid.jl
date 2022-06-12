export Grid, 
    make_grid, 
    basis,
    basis_matrix,
    domain,
    domain_matrix,
    set_domain,
    center, 
    HomeCell, 
    ReciprocalLattice, 
    BrillouinZone, 
    RealLattice, 
    dual_grid,
    transform_grid, 
    snap,
    x_min, x_max, y_min, y_max, z_min, z_max, mins, maxes, 
    HomeCell3D, ReciprocalLattice3D, BrillouinZone3D, RealLattice3D, 
    n_dims, 
    array, 
    invert_grid, 
    shrink, 
    linear_index

"""
WTP provides *non-orthogonal 3D periodic* grids that support

1. Three Indexing schemes and iteration.
2. Fourier transform.

The goals is to provide a uniform interface for dealing with

1. The reciprocal lattice.
2. The Brillouin zone.
3. The crystal lattice.
4. The homecell.

A grid is made of a set of basis vectors and the domain in units of the basis
vectors.  If you subtype this, the concrete type should have these two fields
for things to work out of the box.

```julia
basis::NTuple{3, Vector3}
domain::NTuple{3, NTuple{2, Integer}}
```
"""
abstract type Grid end

"""
    make_grid(T, basis, domain)

# Exmaple: 

Make a grid of type `T` from `basis` and `domain`.

To make a home cell with basis vectors

``a_1 = \\begin{pmatrix}
0\\\\ \\sqrt{3}/2 \\\\ 0
\\end{pmatrix} \\quad
a_2 = \\begin{pmatrix}
1\\\\ 1/2 \\\\ \\sqrt{3}/2 
\\end{pmatrix} \\quad
a_3 = \\begin{pmatrix}
0\\\\ 0 \\\\ 1/2
\\end{pmatrix}``

and domain ``((-2, 1), (-3, 2), (-4, 3))``

```jldoctest grid
julia> homecell_basis = (Vector3(0, sqrt(3)/2, 0), Vector3(1, 1/2, sqrt(3)/2), Vector3(0, 0, 1/2));

julia> homecell = make_grid(HomeCell3D, homecell_basis, ((-2, 1), (-3, 2), (-4, 3)))
type: HomeCell3D
domain: ((-2, 1), (-3, 2), (-4, 3))
basis:
    ket: 0.000, 0.866, 0.000
    ket: 1.000, 0.500, 0.866
    ket: 0.000, 0.000, 0.500

```
"""
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

"""
    basis(g)

Basis vectors for the grid. Depending on the type of grid,
the basis is that of the grid such as real/reciprocal lattice vectors.

A tuple of vectors will be returned, and the vectors returned are always ket
vectors.

Example: 

```jldoctest grid
julia> basis(homecell)

(ket: 0.000, 0.866, 0.000, ket: 1.000, 0.500, 0.866, ket: 0.000, 0.000, 0.500)
```

If it is a matrix that you want, use `basis_matrix` instead for better
performance.
"""
basis(grid::Grid) = grid.basis

"""
    basis_matrix(g)

The basis vector as a matrix, where each column is a basis vector.
The matrix is cached because this has an impact on performance.

Example:

```jldoctest grid
julia> basis_matrix(homecell)

3×3 SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
 0.0       1.0       0.0
 0.866025  0.5       0.0
 0.0       0.866025  0.5
```
"""
basis_matrix(grid::Grid) = grid._cache._basis_matrix


"""
    domain(grid)

The domain of the grid are the extremal miller indices.
The grid includes the end points of the domain.

The result is a tuple of `n_dim` tuples in the form
`((x_min, x_max), (y_min, y_max), (z_min, z_max))`.

Example: 

```jldoctest grid
julia> domain(homecell)

((-2, 1), (-3, 2), (-4, 3))
```

"""
domain(grid::Grid) = grid.domain

_convert_domain_to_matrix(domain) = SMatrix{3,2,Int}(vcat([[d[1] d[2]] for d in domain]...))

"""
    domain_matrix(grid)

Get the domain of the grid as a matrix, where each row 
corresponds to a dimension.

```jldoctest grid
julia> domain_matrix(homecell)

3×2 SMatrix{3, 2, Int64, 6} with indices SOneTo(3)×SOneTo(2):
 -2  1
 -3  2
 -4  3
```
"""
domain_matrix(grid::Grid) = grid._cache._domain_matrix

"""
    set_domain(grid, new_domain)

Create a copy of `grid` with its domain set to `new_domain`.

Example: 

```jldoctest grid
julia> odd_grid = set_domain(homecell, ((-0, 0), (-1, 1), (-2, 2)))

type: HomeCell3D
domain: ((0, 0), (-1, 1), (-2, 2))
basis:
    ket: 0.000, 0.866, 0.000
    ket: 1.000, 0.500, 0.866
    ket: 0.000, 0.000, 0.500
```
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

"""
    expand(grid, factors=[2, 2, 2])

Expand a grid by a factor of `factor`. The returned grid will have the 
save basis vector, but it will be larger.

```jldoctest grid
julia> supercell = expand(homecell, [3, 3, 3])

type: HomeCell3D
domain: ((-6, 5), (-9, 8), (-12, 11))
basis:
    ket: 0.000, 0.866, 0.000
    ket: 1.000, 0.500, 0.866
    ket: 0.000, 0.000, 0.500
```
"""
function expand(grid::Grid, factors::Vector{Int}=[2, 2, 2])
    return set_domain(grid, size_to_domain(factors .* size(grid)))
end

x_min(grid::Grid) = mins(grid)[1]
x_max(grid::Grid) = maxes(grid)[1]
y_min(grid::Grid) = mins(grid)[2]
y_max(grid::Grid) = maxes(grid)[2]
z_min(grid::Grid) = mins(grid)[3]
z_max(grid::Grid) = maxes(grid)[3]


"""
    mins(grid::Grid)

The minimum indices of the grid along each direction as a tuple of integers.

Example: 

```jldoctest grid
julia> mins(homecell)

(-2, -3, -4)
```
"""
mins(grid::Grid) = grid._cache._mins # domain_matrix(grid)[:, 1]

"""
    maxes(grid)

Similar to `mins(grid)`, but gives the maximum indices instead.

Example:

```jldoctest grid
julia> maxes(homecell)

(1, 2, 3)
```

"""
maxes(grid::Grid) = grid._cache._maxes # domain_matrix(grid)[:, 2]

"""
    center(grid)

The center of a grid. This should be a tuple of zeros for a grid that complies 
to the centering convention.  The result will be a tuple of `n_dims` integers.
A tuple instead of a vector is used for performance reasons. 


```jldoctest grid
julia> center(homecell)

(0, 0, 0)
```
"""
center(grid::Grid)::Tuple = grid._cache._center
_compute_center(size::Tuple, domain::Tuple)::Tuple = map(_center_1d, size, domain)
_center_1d(s, d)::Int = iseven(s) ? (sum(d) + 1) ÷ 2 : sum(d) ÷ 2

"""
    n_dims(T)

Gives the number of dimensions of the *type* of grid.

Example:

```jldoctest grid
julia> n_dims(HomeCell3D)

3
```
"""
n_dims(::Type{<:Grid}) = 3

# function center(grid::T) where T <: Grid
#     result = zeros(Int, n_dims(T))
#     s, d = size(grid), domain(grid)
#     for i = 1:n_dims(T)
#         u, l = d[i]
#         result[i] = iseven(s[i]) ? (u + l + 1) ÷ 2 : (u + l) ÷ 2
#     end
#     return result
# end



_compute_size(domain) = Tuple(d[2] - d[1] + 1 for d in domain)

"""
    size(g)

The size of the grid. The result is a tuple of `n_dims(g)` integers.


```jldoctest grid
julia> size(homecell)

(4, 6, 8)
```
"""
Base.size(g::Grid)::Tuple = g._cache._size

"""
    length(g)

The number of grid points in the grid.


```jldoctest grid
julia> length(homecell)

192
```
"""
Base.length(g::Grid)::Int = prod(size(g))

"""
    snap(grid, point)

Snap a cartesian coordinate to a grid point.  

```jldoctest grid
julia> snap(homecell, [0, sqrt(3), 1.49])

GridVector{HomeCell3D}:
    coefficients: [2, 0, 3]
```
"""
snap(grid::Grid, point::AbstractVector)::AbstractGridVector =
    make_grid_vector(grid, (i -> round(Int, i)).((basis_matrix(grid) \ point)))


const Range = AbstractVector{<:Integer}
const NotColon = Union{Range,Integer}

"""
    g[i, j, k]

Indexing a grid with a miller index gives the grid vector corresponding to that
index.

```jldoctest grid
julia> origin = homecell[0, 0, 0]

GridVector{HomeCell3D}:
    coefficients: [0, 0, 0]
```
"""
Base.getindex(g::Grid, i::Integer, j::Integer, k::Integer) = make_grid_vector(g, SVector(i, j, k))

Base.getindex(g::Grid, I::Range, J::Range, K::Range) = [g[i, j, k] for i in I, j in J, k in K]

Base.getindex(g::Grid, i::Integer, J::Range, K::Range) = [g[i, j, k] for j in J, k in K]

"""
    g[x_1:x_2, y, z_1:z_2]

Ranges are supported, and the result will be a multi-dimensional array.

```jldoctest grid
julia> homecell[-1:0, 0, -1:0]

2×2 Matrix{GridVector{HomeCell3D}}:
 GridVector{HomeCell3D}:
    coefficients: [-1, 0, -1]
  GridVector{HomeCell3D}:
    coefficients: [-1, 0, 0]

 GridVector{HomeCell3D}:
    coefficients: [0, 0, -1]
   GridVector{HomeCell3D}:
    coefficients: [0, 0, 0]
```
"""
Base.getindex(g::Grid, I::Range, j::Integer, K::Range) = [g[i, j, k] for i in I, k in K]
Base.getindex(g::Grid, I::Range, J::Range, k::Integer) = [g[i, j, k] for i in I, j in J]

Base.getindex(g::Grid, I::Range, j::Integer, k::Integer) = [g[i, j, k] for i in I]
Base.getindex(g::Grid, i::Integer, J::Range, k::Integer) = [g[i, j, k] for j in J]
Base.getindex(g::Grid, i::Integer, j::Integer, K::Range) = [g[i, j, k] for k in K]

Base.getindex(g::Grid, ::Colon, J::Any, K::Any) = g[x_min(g):x_max(g), J, K]
Base.getindex(g::Grid, I::NotColon, ::Colon, K::Any) = g[I, y_min(g):y_max(g), K]
Base.getindex(g::Grid, I::NotColon, J::NotColon, ::Colon) = g[I, J, z_min(g):z_max(g)]

"""
    linear_index(g::Grid, i::Integer)

Alternative syntax: `g(i)`.

This provides a mapping between 1D indices (which are to be used for matrix
algorithms) and grid vectors.

Example: 

```jldoctest grid
julia> homecell(4)

GridVector{HomeCell3D}:
    coefficients: [-1, 0, 0]

julia> homecell(1:3)

3-element Vector{GridVector{HomeCell3D}}:
 GridVector{HomeCell3D}:
    coefficients: [0, 0, 0]

 GridVector{HomeCell3D}:
    coefficients: [1, 0, 0]

 GridVector{HomeCell3D}:
    coefficients: [-2, 0, 0]
```

The danger of using a one dimensional index is that it can be disassociated with
the object that it's referring to. Therefore, a two way mapping must be
provided. 

```jldoctest grid
julia> linear_index(homecell(4))

4
```
"""
linear_index(g::Grid, i::Integer) = g(i)
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
    iterate(grid)

All grids are iterable. The iteration order is the same order as the one
dimensional index.

```jldoctest grid
julia> collect(homecell)[5] == homecell(5)

true
```

# Example 

For loops

```jldoctest grid
julia> for r in homecell
           r == homecell(1) && println(r)
       end
    
GridVector{HomeCell3D}:
    coefficients: [0, 0, 0]
```

Higher order functions

```jldoctest grid
julia> sum(r->r'*r, homecell)

1779.138438763306
```

Mapping over a grid gives you an `OnGrid` object, which can 
be indexed by grid vectors.

```jldoctest grid
julia> r2 = map(r->norm(r)^2, homecell);
julia> typeof(r2)

WTP.SimpleFunctionOnGrid{HomeCell3D}

julia> r2[homecell[0, 0, 0]]

0
```

One can construct model orbitals this way.

```jldoctest grid
julia> gaussian = map(r->exp(-r'*r), homecell);
julia> contourplot(gaussian[homecell[0, :, :]])
┌────────────────────────────────────────┐  ⠀⠀⠀⠀  
6 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  ┌──┐ 1
  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  │▄▄│  
  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  │▄▄│  
  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  │▄▄│  
  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  │▄▄│  
  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣀⣀⣀⣀⠤⠤⣀⡀⠀⠀⠀│  │▄▄│  
  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⣀⡠⠤⠤⠔⠒⠒⠒⠉⠉⠉⠉⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀│  │▄▄│  
  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⡠⠤⠒⠉⢀⣀⡠⠤⠔⠒⠒⠒⠒⠒⠒⠒⠒⢄⠀⠀⠀⠀⠀⠀⢱⠀⠀⠀│  │▄▄│  
  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⡠⠊⠁⠀⢀⠤⠒⠁⣀⡠⠤⠒⠒⠉⠑⠒⠢⠤⠤⣀⠀⠉⠢⡀⠀⠀⠀⢸⠀⠀⠀│  │▄▄│  
  │⠀⠀⠀⠀⠀⠀⠀⠀⠸⡀⠀⠀⠀⠉⠢⣀⠈⠑⠒⠒⠤⠤⣀⣀⠤⠤⠒⠊⠉⡠⠔⠊⠁⠀⣀⠔⠉⠀⠀⠀│  │▄▄│  
  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠑⠤⠤⠤⠤⠤⠤⠤⠤⠔⠒⠊⠉⢀⡠⠔⠒⠉⠀⠀⠀⠀⠀⠀│  │▄▄│  
  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⣀⣀⣀⡠⠤⠤⠤⠔⠒⠒⠊⠉⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  │▄▄│  
  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠑⠒⠒⠉⠉⠉⠉⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  │▄▄│  
  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  │▄▄│  
1 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│  └──┘ 0
  └────────────────────────────────────────┘  ⠀⠀⠀⠀  
  ⠀1⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀8⠀  ⠀⠀⠀⠀  
```

The result does not look isotropic even though our Gaussian is. This is because
the grid is not orthogonal. Keep this in mind for all the pictures.
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

If we apply a Fourier transform on a function defined on a grid, the result is a
function defined on a different grid, which we refer to as the dual grid.

```jldoctest grid
julia> guassian_reciprocal = fft(gaussian);
julia> reciprocal_lattice = grid(guassian_reciprocal)

type: ReciprocalLattice3D
domain: ((-2, 1), (-3, 2), (-4, 3))
basis:
    ket: -0.907, 1.814, -0.000
    ket: 1.047, 0.000, -0.000
    ket: -1.360, 0.000, 1.571
```

Example: 

```jldoctest grid
julia> reciprocal_lattice = transform_grid(homecell)

type: ReciprocalLattice3D
domain: ((-2, 1), (-3, 2), (-4, 3))
basis:
    ket: -0.907, 1.814, -0.000
    ket: 1.047, 0.000, -0.000
    ket: -1.360, 0.000, 1.571
```

The resulting basis vectors should satisfy aᵢᵀ bᵢ = 2π / sᵢ, where s is the size
of the grid (number of grid points in each direction).

```jldoctest grid
julia> basis_matrix(reciprocal_lattice)' * basis_matrix(homecell) * diagm([size(homecell)...])
3×3 Matrix{Float64}:
 6.28319   0.0          0.0
 0.0       6.28319      0.0
 0.0      -6.31568e-17  6.28319
```

`HomeCell3D` and `ReciprocalLattice3D` are dual grids and
`BrillouinZone3D` and `RealLattice3D` are dual grids.

"""
transform_grid(grid::T) where {T<:Grid} =
    let A = basis_matrix(grid) * diagm([size(grid)...]), B = 2 * pi * inv(A)'
        make_grid(dual_grid(T), matrix_to_vector3(B), size_to_domain(size(grid)))
    end


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


dual_grid(::Type{HomeCell3D}) = ReciprocalLattice3D
dual_grid(::Type{ReciprocalLattice3D}) = HomeCell3D
dual_grid(::Type{BrillouinZone3D}) = RealLattice3D
dual_grid(::Type{RealLattice3D}) = BrillouinZone3D

# inverse_grid(::Type{HomeCell3D}) = RealLattice3D
# inverse_grid(::Type{RealLattice3D}) = HomeCell3D
# inverse_grid(::Type{BrillouinZone3D}) = ReciprocalLattice3D
# inverse_grid(::Type{ReciprocalLattice3D}) = BrillouinZone3D

Base.:(==)(g_1::T, g_2::T) where {T<:Grid} = basis(g_1) == basis(g_2) && domain(g_1) == domain(g_2)