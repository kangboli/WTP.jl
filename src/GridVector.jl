export AbstractGridVector,
    GridVector,
    wrapped,
    overflow,
    reset_overflow,
    has_overflow,
    grid_vector_constructor,
    add_overflow,
    cartesian,
    linear_index

abstract type AbstractGridVector{G<:Grid} <: LinearCombination end

grid(grid_vector::AbstractGridVector) = grid_vector.grid
grid!(grid_vector::AbstractGridVector, new_grid) = grid_vector.grid = new_grid

"""
The grid vector modulo the domain of the grid.

For grid vectors within the grid domain, this just gives the _coefficients.
For those outside the grid domain, this move the vector by multiples of 
the grid domain to bring it into the grid domain.
"""
wrapped(grid_vector::AbstractGridVector) = [
    mod((c - l), r - l + 1) + l for
    (c, (l, r)) in zip(coefficients(grid_vector), domain(grid(grid_vector)))
]

"""
This gives the number of times we have to translate a vector by a grid domain to
bring them into the grid domain.
"""
overflow(grid_vector::AbstractGridVector) = [
    fld((c - l), r - l + 1) for
    (c, (l, r)) in zip(coefficients(grid_vector), domain(grid(grid_vector)))
]


"""
The gives a image of the grid_vector within the grid domain.
"""
function reset_overflow(grid_vector::AbstractGridVector)
    @set grid_vector._coefficients = wrapped(grid_vector)
end

"""
This moves the grid_vector by some multiples of the grid domain.
"""
function add_overflow(grid_vector::AbstractGridVector, overflow)
    @set grid_vector._coefficients = [
        c + o * s for
        (o, c, s) in zip(overflow, coefficients(grid_vector), size(grid(grid_vector)))
    ]
end

"""
Whether there is any overflow.
"""
has_overflow(grid_vector::AbstractGridVector) =
    let (fx, fy, fz) = overflow(grid_vector)
        fx != 0 || fy != 0 || fz != 0
    end


struct GridVector{T} <: AbstractGridVector{T}
    grid::T
    _coefficients::AbstractVector
    ket::Bool
end

"""
Return the basis as they are internally stored (as kets).
"""
_basis(grid_vector::GridVector) = basis(grid(grid_vector))

function add(grid_vec_1::GridVector, grid_vec_2::GridVector)
    ket(grid_vec_1) == ket(grid_vec_2) || error("Adding a bra to a ket.")
    grid(grid_vec_1) == grid(grid_vec_2) || error("Grid vectors defined on different grid.")
    return typeof(grid_vec_1)(
        grid(grid_vec_1),
        coefficients(grid_vec_1) + coefficients(grid_vec_2),
        ket(grid_vec_1),
    )
end

negate(grid_vec_1::GridVector) =
    typeof(grid_vec_1)(grid(grid_vec_1), -coefficients(grid_vec_1), ket(grid_vec_1))
mul(s::Number, l1::GridVector) = typeof(l1)(grid(l1), s * coefficients(l1), ket(l1))
minus(grid_vec_1::GridVector, grid_vec_2::GridVector) = add(grid_vec_1, negate(grid_vec_2))

grid_vector_constructor(g::T, _coefficients) where {T<:Grid} =
    GridVector{T}(g, SVector(_coefficients...), true)

function Base.show(io::IO, grid_vector::GridVector)
    print(io, "grid: $(typeof(grid(grid_vector)))\n")
    print(io, "_coefficients: $(coefficients(grid_vector))\n")
end

Base.:(==)(grid_vec_1::GridVector, grid_vec_2::GridVector) =
    coefficients(grid_vec_1) == coefficients(grid_vec_2) &&
    ## TODO: Put this back and figure out if it breaks things.
    # grid(grid_vec_1) == grid(grid_vec_2) &&
    ket(grid_vec_1) == ket(grid_vec_2)

Base.hash(grid_vec::GridVector) = hash(coefficients(grid_vec)) + hash(ket(grid_vec))

"""
    cartesian(grid_vec)

The Cartesian coordinates of a grid vector.
"""
cartesian(grid_vec::GridVector)::Vector{Number} = basis_transform(
    coefficients(grid_vec), basis(grid_vec), CARTESIAN_BASIS)

"""
Convert a grid vector to an integer as a linear index.
Indexing the underlying grid with this integer gives back
the grid vector.
"""
linear_index(grid_vec::GridVector) = let sizes = size(grid(grid_vec))
    three_to_one(miller_to_standard(sizes, coefficients(grid_vec), [0, 0, 0])..., sizes)
end