export AbstractGridVector,
    GridVector,
    grid,
    set_grid,
    wrapped,
    overflow,
    reset_overflow,
    has_overflow,
    make_grid_vector,
    add_overflow,
    coordinates,
    linear_index,
    braket_

"""
A grid vector is a point on our grid. It is primarily used to index functions
defined on a grid, but it also behaves like a grid point, which means that you
can do arithmatics with them. 
"""
abstract type AbstractGridVector{G<:Grid} <: LinearCombination end

struct GridVector{T <: Grid} <: AbstractGridVector{T}
    grid::T
    _coefficients::SVector{3, Int}
    ket::Bool
end

"""
    make_grid_vector(g, coefficients)

Make a grid vector that is at g[coefficients...].

```jldoctest grid_vector
julia> reciprocal_basis = CARTESIAN_BASIS;

julia> sizes = (4, 4, 4);

julia> lattice = make_grid(ReciprocalLattice3D, reciprocal_basis, size_to_domain(sizes));

julia> v = make_grid_vector(lattice, [0, 0, 0])
GridVector{ReciprocalLattice3D}:
    coefficients: [0, 0, 0]
```

This is what is behind indexing a grid

```jldoctest grid_vector
julia> v = lattice[0, 0, 0]

GridVector{ReciprocalLattice3D}:
    coefficients: [0, 0, 0]
```
"""
make_grid_vector(g::T, _coefficients::AbstractVector) where {T<:Grid} =
    GridVector{T}(g, _coefficients, true)

"""
    grid(grid_vector)

Each grid vector is defined on a grid.
This gets the grid on which the grid vector is defined.

Example:

```jldoctest grid_vector
julia> grid(v)
type: ReciprocalLattice3D
domain: ((-2, 1), (-2, 1), (-2, 1))
basis:
    ket: 1.000, 0.000, 0.000
    ket: 0.000, 1.000, 0.000
    ket: 0.000, 0.000, 1.000
```
"""
grid(grid_vector::AbstractGridVector) = grid_vector.grid


"""
    set_grid(grid_vector, new_grid)

Put the grid_vector on the `new_grid` while preserving the coefficients.

Example:

```@jldoctest grid_vector
julia> ṽ = set_grid(v, expand(lattice));

julia> grid(ṽ)

type: ReciprocalLattice3D
domain: ((-4, 3), (-4, 3), (-4, 3))
basis:
    ket: 1.000, 0.000, 0.000
    ket: 0.000, 1.000, 0.000
    ket: 0.000, 0.000, 1.000
```
"""
set_grid(grid_vector::AbstractGridVector, new_grid) = @set grid_vector.grid = new_grid

miller_to_standard(grid_vector::AbstractGridVector, offsets::Tuple) =
    miller_to_standard(size(grid(grid_vector)), tuple(coefficients(grid_vector)...), offsets)

# miller_to_standard(grid_vector::AbstractGridVector{T}, offsets::AbstractGridVector{T}) where {T} =
#     miller_to_standard(grid_vector, coefficients(offsets))



"""
    overflow(grid_vector)

A grid vector is allowed to be outside its underlying grid. One can 
move the grid vector by integer multiples of the grid domain so that 
they reside in the grid.

The overflow gives the number of times we have to translate a vector by a grid
domain to bring them into the grid domain.

Example:

```@jldoctest grid_vector
julia> overflow(lattice[0, 0, 0])
3-element Vector{Int64}:
 0
 0
 0

julia> overflow(lattice[3, -3, -3])
3-element Vector{Int64}:
  1
 -1
 -1
```
"""
function overflow(grid_vector::AbstractGridVector) 
    g = grid(grid_vector)
    D = n_dims(typeof(g))
    result = zeros(Int, D)
    for i = 1:D
        result[i] = overflow_1d(coefficients(grid_vector)[i], mins(g)[i], maxes(g)[i])
    end
    return result
end
overflow_1d(c, l, r) = fld((c - l), r - l + 1)

"""
    has_overflow(grid_vector)

Whether there is any overflow.

Example:

```@jldoctest grid_vector
julia> has_overflow(lattice[0, 0, 5])
true

julia> has_overflow(lattice[0, 0, 0])
false
```
"""
has_overflow(grid_vector::AbstractGridVector) =
    any(f -> f != 0, overflow(grid_vector))

"""
    wrapped(grid_vector)

The grid vector modulo the domain of the grid. One can think of this a 
the remainder of the grid vector divided by the grid.

For grid vectors within the grid domain, this just gives the _coefficients.  For
those outside the grid domain, this move the vector by multiples of the grid
domain to bring it into the grid domain. 

Example:

```@jldoctest grid_vector
julia> wrapped(lattice[0, 1, 0])
3-element Vector{Int64}:
 0
 1
 0

julia> wrapped(lattice[3, -3, -3])
3-element Vector{Int64}:
 -1
  1
  1
```
"""
function wrapped(grid_vector::AbstractGridVector) 
    g = grid(grid_vector)
    D = n_dims(typeof(g))
    result = zeros(Int, D)
    for i = 1:D
        result[i] = wrapped_1d(coefficients(grid_vector)[i], mins(g)[i], maxes(g)[i])
    end
    return result
end
wrapped_1d(c, l, r) = mod((c - l), r - l + 1) + l 

function wrapped(grid::Grid, point)
    snapped = snap(grid, point)
    residual = point - coordinates(snapped)
    return residual + coordinates(reset_overflow(snapped))
end

# An alternative implementation.
# function wrapped(grid_vector::AbstractGridVector)
#     c = coefficients(grid_vector)
#     d = domain(grid(grid_vector))
#     result = zeros(Int, 3)
#     for i = 1:3
#         (l, r) = d[i]
#         result[i] = mod((c[i] - l), r - l + 1) + l
#     end
#     return result
# end
"""
    reset_overflow(grid_vector)

The gives a image of the grid_vector within the grid domain.

Example:

```@jldoctest grid_vector
julia> reset_overflow(lattice[0, 0, 5])
GridVector{ReciprocalLattice3D}:
    coefficients: [0, 0, 1]
```
"""
function reset_overflow(grid_vector::AbstractGridVector) 
    make_grid_vector(grid(grid_vector), wrapped(grid_vector))
end

"""
    add_overflow(grid_vector, overflow)

This moves the grid_vector by some multiples of the grid domain.

Example:

```@jldoctest grid_vector
julia> add_overflow(lattice[0, 0, 1], [0, 0, 1])
GridVector{ReciprocalLattice3D}:
    coefficients: [0, 0, 5]
```
"""
function add_overflow(grid_vector::AbstractGridVector, overflow) 
    new_coefficients = map((o, c, s)->c + o * s, overflow, coefficients(grid_vector), SVector(size(grid(grid_vector))))
    make_grid_vector(grid(grid_vector), new_coefficients)
end


"""
Return the basis as they are internally stored (as kets).
"""
_basis(grid_vector::AbstractGridVector) = basis(grid(grid_vector))

"""
    add(grid_vector_1, grid_vector_2)

Add two grid vectors. Can also write `grid_vector_1 + grid_vector_2`.

Example:

```@jldoctest grid_vector
julia> lattice[0, 0, 1] + lattice[-1, 0, 0]
GridVector{ReciprocalLattice3D}:
    coefficients: [-1, 0, 1]
```
"""
function add(grid_vector_1::T, grid_vector_2::T) where T <: AbstractGridVector
    # ket(grid_vector_1) == ket(grid_vector_2) || error("Adding a bra to a ket.")
    grid(grid_vector_1) == grid(grid_vector_2) || error("Grid vectors defined on different grid.")
    new_coefficients = coefficients(grid_vector_1) + coefficients(grid_vector_2)
    return T(grid(grid_vector_1), new_coefficients, ket(grid_vector_1))
end

"""
    negate(grid_vector_1)

Can also write `-grid_vector_1`.

Example:

```@jldoctest grid_vector
julia> -lattice[0, 0, 1]
GridVector{ReciprocalLattice3D}:
    coefficients: [0, 0, -1]
```
"""
negate(grid_vector_1::T) where T <: AbstractGridVector =
    T(grid(grid_vector_1), -coefficients(grid_vector_1), ket(grid_vector_1))

"""
    mul(s, grid_vector)

Scale a grid vector. Can also write `s * grid_vector` or `grid_vector * s`.
"""
mul(s::Int, l1::T) where T <: AbstractGridVector = T(grid(l1), s * coefficients(l1), ket(l1))
# minus(grid_vector_1::T, grid_vector_2::T) where T <: AbstractGridVector = add(grid_vector_1, negate(grid_vector_2))

function Base.show(io::IO, grid_vector::AbstractGridVector)
    print(io, "$(typeof(grid_vector)):\n")
    # print(io, "grid: $(typeof(grid(grid_vector)))\n" |> indent)
    print(io, indent("coefficients: $(coefficients(grid_vector))") * "\n")
end

Base.:(==)(grid_vector_1::T, grid_vector_2::T) where T <: AbstractGridVector =
    coefficients(grid_vector_1) == coefficients(grid_vector_2) # && grid(grid_vector_1) == grid(grid_vector_2) 

Base.hash(grid_vector::AbstractGridVector)::UInt = linear_index(grid_vector) |> abs

"""
    coordinates(grid_vector)

The Cartesian coordinates of a grid vector.

Example:

```@jldoctest grid_vector
julia> coordinates(lattice[1, 0, 1])
3-element Vector{Number}:
 1.0
 0.0
 1.0
```
"""
coordinates(grid_vector::AbstractGridVector)::Vector{Number} = let g = grid(grid_vector)
    basis_matrix(g) * (coefficients(grid_vector) + SVector(center(g)))
end
# cartesian(grid_vec::AbstractGridVector)::Vector{Number} = basis_transform(
#     coefficients(grid_vec), basis(grid_vec), CARTESIAN_BASIS)

"""
    norm(grid_vector)

The 2-norm of a grid vector

Example:

```@jldoctest grid_vector
julia> norm(lattice[1, 1, 1])
1.7320508075688772
```
"""
LinearAlgebra.norm(grid_vector::AbstractGridVector) = norm(coordinates(grid_vector))

"""
    linear_index(grid_vector)


Convert a grid vector to an integer as a linear index.
Indexing the underlying grid with this integer gives back
the grid vector.
"""
linear_index(grid_vector::AbstractGridVector) =
    let g = grid(grid_vector), s = size(g)
        three_to_one(miller_to_standard(s, tuple(coefficients(grid_vector)...), center(g))..., s)
    end

"""
    braket(grid_vector_1, grid_vector_2) 

Compute the inner product of two grid vectors. Can also write
`grid_vector_1' * grid_vector_2`

PS: a `braket_` function is created as a dummy anchor for documentation.

Example:

```@jldoctest grid_vector
julia> homecell = transform_grid(lattice);

julia> homecell[0, 0, 1]' * lattice[1, 0, 1]
1.5707963267948966
```
"""
braket_(grid_vector_1::AbstractGridVector, grid_vector_2::AbstractGridVector) = braket(grid_vector_1, grid_vector_2)
