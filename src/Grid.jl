export Grid, domain, domain!, basis, dual_grid, 
HomeCell, ReciprocalLattice, BrillouinZone, RealLattice, transform_grid, snap

"""
Abstract Grid. 

This abstract type provides indexing with Miller indices (slightly generalized)
and iteration.
"""

abstract type Grid end

"""
The domain of the grid are the extremal miller indices.
The grid includes the end points of the domain.
"""
domain(grid::Grid) = grid.domain
domain!(grid::Grid, new_domain) = grid.domain = new_domain

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

"""
    g[0, 0, 1]

Indexing a grid with a miller index gives the grid vector corresponding to that
index.
"""
Base.getindex(g::Grid, i::Integer, j::Integer, k::Integer) = grid_vector_constructor(g, [i, j, k])

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
Base.getindex(g::Grid, linear_indices::AbstractVector{<:Integer}) = [g[i] for i in linear_indices]

"""
Iterate over the grid gives a sequence of grid vectros that goes through each grid point.
"""
function Base.iterate(grid::Grid)
    n_x, n_y, n_z = size(grid)
    (n_x < 1 || n_y < 1 || n_z < 1) && return nothing
    miller = Iterators.product(
        -div(n_x, 2)+1:div(n_x, 2),
        -div(n_y, 2)+1:div(n_y, 2),
        -div(n_z, 2)+1:div(n_z, 2))
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

struct HomeCell <: Grid
    basis::Tuple{Vector3, Vector3, Vector3}
    domain::Tuple{Tuple{Integer, Integer}, Tuple{Integer, Integer}, Tuple{Integer, Integer}}
end

struct ReciprocalLattice <: Grid
    basis::Tuple{Vector3, Vector3, Vector3}
    domain::Tuple{Tuple{Integer, Integer}, Tuple{Integer, Integer}, Tuple{Integer, Integer}}
end

struct BrillouinZone <: Grid
    basis::Tuple{Vector3, Vector3, Vector3}
    domain::Tuple{Tuple{Integer, Integer}, Tuple{Integer, Integer}, Tuple{Integer, Integer}}
end

struct RealLattice <: Grid
    basis::Tuple{Vector3, Vector3, Vector3}
    domain::Tuple{Tuple{Integer, Integer}, Tuple{Integer, Integer}, Tuple{Integer, Integer}}
end

dual_grid(::Type{HomeCell}) = ReciprocalLattice
dual_grid(::Type{ReciprocalLattice}) =  HomeCell
dual_grid(::Type{BrillouinZone}) = RealLattice
dual_grid(::Type{RealLattice}) = BrillouinZone