export Grid, domain, domain!, basis, dual_grid, 
HomeCell, ReciprocalLattice, BrillouinZone, RealLattice, transform_grid, snap,
x_min, x_max, y_min, y_max, z_min, z_max, mins, maxes,
HomeCell3D, ReciprocalLattice3D, BrillouinZone3D, RealLattice3D

"""
Abstract Grid. 

This abstract type provides indexing with Miller indices (slightly generalized)
and iteration.

Centering convention:

    The convention for grid centered at the origin:
        - For an even grid: the domain is -N, ..., 0, ..., N-1.
        - For an odd grid: the domain is -N, ..., 0, ..., N.

This convention matches the input files of QE for the Brillouin zone (but not
its internals and outputs). This also matches the reciprocal lattice 
and possibly the homecell in QE.

If a grid does not comply to the centering convention, it is considered a
translated grid.
"""

abstract type Grid end

"""
The domain of the grid are the extremal miller indices.
The grid includes the end points of the domain.
"""
domain(grid::Grid) = grid.domain
domain!(grid::Grid, new_domain) = grid.domain = new_domain
x_min(grid::Grid) = domain(grid)[1][1]
x_max(grid::Grid) = domain(grid)[1][2]
y_min(grid::Grid) = domain(grid)[2][1]
y_max(grid::Grid) = domain(grid)[2][2]
z_min(grid::Grid) = domain(grid)[3][1]
z_max(grid::Grid) = domain(grid)[3][2]

mins(grid::Grid) = [d[1] for d in domain(grid)]
maxes(grid::Grid) = [d[2] for d in domain(grid)]
array(grid::Grid) = grid[:,:,:]

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
Base.getindex(g::Grid, i::Integer) = grid_vector_constructor(g,
    standard_to_miller(size(g), one_to_three(i, size(g)), [0,0,0]))

"""
Indexing a grid with a list of linear indices gives a list of grid vector.
"""
Base.getindex(g::Grid, linear_indices::Range) = [g[i] for i in linear_indices]

"""
Iterate over the grid gives a sequence of grid vectros that goes through each grid point.
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
FFT of the grid.

The resulting basis vectors should satisfy
aᵢ bᵢ = 2π / sᵢ,
where s is the size of the grid.
"""
transform_grid(grid::T) where T <: Grid = 
    let A = vector3_to_matrix(basis(grid)) * diagm([size(grid)...]), B = 2 * pi * inv(A)'
        dual_grid(T)(matrix_to_vector3(B), size_to_domain(size(grid)))
    end

"""
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

function Base.show(io::IO, grid::Grid)
    println(io, string(grid))
end


"""
There are four basic grids (two pairs of dual grids) in Condensed Phase on which
most physical concepts are defined. Generic grid operations should be defined
on the abstract type `Grid`. Operations specific to each space should be defined
on the concrete grid.

This approach is to be contrasted with allocating a Fortran array and program  
functions as array manipulation.
"""

abstract type HomeCell <: Grid end
struct HomeCell3D <: HomeCell
    basis::Tuple{Vector3, Vector3, Vector3}
    domain::Tuple{Tuple{Integer, Integer}, Tuple{Integer, Integer}, Tuple{Integer, Integer}}
end

abstract type ReciprocalLattice <: Grid end 
struct ReciprocalLattice3D <: ReciprocalLattice
    basis::Tuple{Vector3, Vector3, Vector3}
    domain::Tuple{Tuple{Integer, Integer}, Tuple{Integer, Integer}, Tuple{Integer, Integer}}
end

abstract type BrillouinZone <: Grid end
struct BrillouinZone3D <: BrillouinZone
    basis::Tuple{Vector3, Vector3, Vector3}
    domain::Tuple{Tuple{Integer, Integer}, Tuple{Integer, Integer}, Tuple{Integer, Integer}}
end

abstract type RealLattice <: Grid end
struct RealLattice3D <: RealLattice
    basis::Tuple{Vector3, Vector3, Vector3}
    domain::Tuple{Tuple{Integer, Integer}, Tuple{Integer, Integer}, Tuple{Integer, Integer}}
end

dual_grid(::Type{HomeCell3D}) = ReciprocalLattice3D
dual_grid(::Type{ReciprocalLattice3D}) =  HomeCell3D
dual_grid(::Type{BrillouinZone3D}) = RealLattice3D
dual_grid(::Type{RealLattice3D}) = BrillouinZone3D