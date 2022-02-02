export AbstractGridVector,
    GridVector,
    wrapped,
    overflow,
    reset_overflow,
    has_overflow,
    grid_vector_constructor,
    add_overflow,
    cartesian,
    linear_index,
    grid

abstract type AbstractGridVector{G<:Grid} <: LinearCombination end

grid(grid_vector::AbstractGridVector) = grid_vector.grid
grid!(grid_vector::AbstractGridVector, new_grid) = grid_vector.grid = new_grid

"""
    wrapped(grid_vector)

The grid vector modulo the domain of the grid.

For grid vectors within the grid domain, this just gives the _coefficients.
For those outside the grid domain, this move the vector by multiples of 
the grid domain to bring it into the grid domain.
"""
# wrapped(grid_vector::AbstractGridVector) = 
#     coefficients(grid_vector) + overflow(grid_vector) .* [size(grid(grid_vector))...]

function wrapped(grid_vector::AbstractGridVector) 
    c = coefficients(grid_vector)
    d = domain(grid(grid_vector))
    result = zeros(Int, 3)
    for i = 1:3
        (l, r) = d[i]
        result[i] = mod((c[i] - l), r - l + 1) + l
    end
    return result
end

"""
This gives the number of times we have to translate a vector by a grid domain to
bring them into the grid domain.
"""
function overflow(grid_vector::AbstractGridVector) 
    c = coefficients(grid_vector)
    d = domain(grid(grid_vector))
    result = zeros(Int, 3)
    for i = 1:3
        (l, r) = d[i]
        result[i] = fld((c[i] - l), r - l + 1)
    end
    return result
end


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
    any(f->f!=0, overflow(grid_vector))


struct GridVector{T} <: AbstractGridVector{T}
    grid::T
    _coefficients::AbstractVector
    ket::Bool
end

"""
Return the basis as they are internally stored (as kets).
"""
_basis(grid_vector::AbstractGridVector) = basis(grid(grid_vector))

function add(grid_vector_1::AbstractGridVector, grid_vector_2::AbstractGridVector)
    ket(grid_vector_1) == ket(grid_vector_2) || error("Adding a bra to a ket.")
    grid(grid_vector_1) == grid(grid_vector_2) || error("Grid vectors defined on different grid.")
    return typeof(grid_vector_1)(
        grid(grid_vector_1),
        coefficients(grid_vector_1) + coefficients(grid_vector_2),
        ket(grid_vector_1),
    )
end

negate(grid_vector_1::AbstractGridVector) =
    typeof(grid_vector_1)(grid(grid_vector_1), -coefficients(grid_vector_1), ket(grid_vector_1))
mul(s::Number, l1::AbstractGridVector) = typeof(l1)(grid(l1), s * coefficients(l1), ket(l1))
minus(grid_vector_1::AbstractGridVector, grid_vector_2::AbstractGridVector) = add(grid_vector_1, negate(grid_vector_2))

grid_vector_constructor(g::T, _coefficients) where {T<:Grid} =
    GridVector{T}(g, SVector(_coefficients...), true)

function Base.show(io::IO, grid_vector::AbstractGridVector)
    print(io, "$(typeof(grid_vector)):\n")
    # print(io, "grid: $(typeof(grid(grid_vector)))\n" |> indent)
    print(io, indent("_coefficients: $(coefficients(grid_vector))")*"\n")
end

Base.:(==)(grid_vector_1::AbstractGridVector, grid_vector_2::AbstractGridVector) =
    coefficients(grid_vector_1) == coefficients(grid_vector_2) &&
    ## TODO: Put this back and figure out if it breaks things.
    # grid(grid_vec_1) == grid(grid_vec_2) &&
    ket(grid_vector_1) == ket(grid_vector_2)

Base.hash(grid_vector::AbstractGridVector) = hash(coefficients(grid_vector)) + hash(ket(grid_vector))

"""
    cartesian(grid_vector)

The Cartesian coordinates of a grid vector.
"""
cartesian(grid_vector::AbstractGridVector)::Vector{Number} = 
    vector3_to_matrix(basis(grid_vector)) * coefficients(grid_vector)
# cartesian(grid_vec::AbstractGridVector)::Vector{Number} = basis_transform(
#     coefficients(grid_vec), basis(grid_vec), CARTESIAN_BASIS)

LinearAlgebra.norm(grid_vector::AbstractGridVector) = norm(cartesian(grid_vector))

"""
Convert a grid vector to an integer as a linear index.
Indexing the underlying grid with this integer gives back
the grid vector.
"""
linear_index(grid_vector::AbstractGridVector) = let sizes = size(grid(grid_vector))
    three_to_one(miller_to_standard(sizes, coefficients(grid_vector), [0, 0, 0])..., sizes)
end
