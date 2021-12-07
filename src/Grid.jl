export Grid, domain, domain!, basis, dual_grid, 
HomeCell, ReciprocalLattice, BrillouinZone, RealLattice, transform_grid, snap,
x_min, x_max, y_min, y_max, z_min, z_max, mins, maxes,
HomeCell3D, ReciprocalLattice3D, BrillouinZone3D, RealLattice3D, n_dims, array, invert_grid

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

"""
    domain(grid)

The domain of the grid are the extremal miller indices.
The grid includes the end points of the domain.
"""
domain(grid::Grid) = grid.domain

"""
    domain!(grid, new_grid)

Set the domain of `grid` to `new_grid`.
"""
domain!(grid::Grid, new_domain) = @set grid.domain = new_domain

x_min(grid::Grid) = domain(grid)[1][1]
x_max(grid::Grid) = domain(grid)[1][2]
y_min(grid::Grid) = domain(grid)[2][1]
y_max(grid::Grid) = domain(grid)[2][2]
z_min(grid::Grid) = domain(grid)[3][1]
z_max(grid::Grid) = domain(grid)[3][2]

"""
    mins(grid::Grid)

The minimum indices of the grid along each direction.

Example: 

```julia
homecell = HomeCell3D(CARTESIAN_BASIS, ((-2, 1), (-2, 1), (-2, 1)))
mins(homecell) # gives [-2, -2, -2]
```
"""
mins(grid::Grid) = [d[1] for d in domain(grid)]

"""
    maxes(grid)

Similar to `mins(grid)`, but gives the maximum indices instead.
"""
maxes(grid::Grid) = [d[2] for d in domain(grid)]

"""
    array(grid)

Convert the grid to a multidimensional array.


Example:

```julia
homecell = HomeCell3D(CARTESIAN_BASIS, ((-2, 1), (-2, 1), (-2, 1)))
array(homecell)[1, 1, 1]
# grid: HomeCell3D
# _coefficients: [-2, -2, -2]
```
"""
array(grid::Grid) = grid[:,:,:]

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

The center of a grid. This should be 0 for a grid that complies 
to the centering convention.
"""
center(grid::Grid) = [(u + l + 1) ÷ 2 for (l, u) in domain(grid)]

"""
    basis(g)

Basis vectors for the grid. Depending on the type of grid,
the basis is that of the grid such as real/reciprocal lattice vectors.

The vectors returned are always ket vectors.
"""
basis(grid::Grid) = grid.basis

"""
    dual_grid(T)

The type of the dual grid of T.
"""
dual_grid(::Type{T}) where T = T

"""
    size(g)

The size of the grid. The result is a tuple of integers.
"""
Base.size(g::Grid) = Tuple(d[2]-d[1]+1 for d in domain(g)) 
Base.length(g::Grid) = prod(size(g))

const Range = AbstractVector{<:Integer}
const NotColon = Union{Range, Integer}
"""
    g[0, 0, 1]

Indexing a grid with a miller index gives the grid vector corresponding to that
index.
"""
Base.getindex(g::Grid, i::Integer, j::Integer, k::Integer) = grid_vector_constructor(g, [i, j, k])

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
(g::Grid)(i::Integer) = grid_vector_constructor(g,
    standard_to_miller(size(g), one_to_three(i, size(g)), [0,0,0]))

"""
Indexing a grid with a list of linear indices gives a list of grid vector.
"""
(g::Grid)(linear_indices::Range) = [g(i) for i in linear_indices]
(g::Grid)(::Colon) = [g(i) for i in 1:length(g)]

"""
Iterate over the grid gives a sequence of grid vectors that goes through each grid point.
"""
function Base.iterate(grid::Grid)

    miller = Iterators.product([l:u for (l, u) in zip(mins(grid), maxes(grid))]...)
    first, miller_state = iterate(miller)
    return (grid_vector_constructor(grid, [first...]), (miller, miller_state))
end

function Base.iterate(grid::Grid, state)
    miller, miller_state = state
    next_miller = iterate(miller, miller_state)
    next_miller === nothing && return nothing
    next, miller_state = next_miller
    return (grid_vector_constructor(grid, [next...]), (miller, miller_state))
end

"""
    transform_grid(g)

FFT of the grid.

The resulting basis vectors should satisfy
aᵢᵀ bᵢ = 2π / sᵢ,
where s is the size of the grid (number of grid points in each direction).
"""
transform_grid(grid::T) where T <: Grid = 
    let A = vector3_to_matrix(basis(grid)) * diagm([size(grid)...]), B = 2 * pi * inv(A)'
        dual_grid(T)(matrix_to_vector3(B), size_to_domain(size(grid)))
    end

invert_grid(grid::T) where T <: Grid = 
    inverse_grid(T)(
        Tuple(2b / s for (b, s) in zip(basis(grid), size(grid))), 
        domain(grid)
    )

"""
    snap(grid, point)

Snap a coordinate to a grid point.  
"""
snap(grid::Grid, point::AbstractVector) =
    grid_vector_constructor(grid, Int.(round.(vector3_to_matrix(basis(grid)) \ point)))


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
    basis::Tuple{Vector3, Vector3, Vector3}
    domain::Tuple{Tuple{Integer, Integer}, Tuple{Integer, Integer}, Tuple{Integer, Integer}}
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
    basis::Tuple{Vector3, Vector3, Vector3}
    domain::Tuple{Tuple{Integer, Integer}, Tuple{Integer, Integer}, Tuple{Integer, Integer}}
end

abstract type BrillouinZone <: Grid end
"""
The 3D Brillouin zone. The usage is the same of `HomeCell3D`.
"""
struct BrillouinZone3D <: BrillouinZone
    basis::Tuple{Vector3, Vector3, Vector3}
    domain::Tuple{Tuple{Integer, Integer}, Tuple{Integer, Integer}, Tuple{Integer, Integer}}
end

abstract type RealLattice <: Grid end
"""
The 3D crystal lattice. This is the dual grid of `BrillouinZone3D`.
"""
struct RealLattice3D <: RealLattice
    basis::Tuple{Vector3, Vector3, Vector3}
    domain::Tuple{Tuple{Integer, Integer}, Tuple{Integer, Integer}, Tuple{Integer, Integer}}
end

dual_grid(::Type{HomeCell3D}) = ReciprocalLattice3D
dual_grid(::Type{ReciprocalLattice3D}) = HomeCell3D
dual_grid(::Type{BrillouinZone3D}) = RealLattice3D
dual_grid(::Type{RealLattice3D}) = BrillouinZone3D

inverse_grid(::Type{HomeCell3D}) = RealLattice3D
inverse_grid(::Type{RealLattice3D}) = HomeCell3D
inverse_grid(::Type{BrillouinZone3D}) = ReciprocalLattice3D
inverse_grid(::Type{ReciprocalLattice3D}) = BrillouinZone3D
