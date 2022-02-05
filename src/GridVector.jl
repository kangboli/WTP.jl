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
    cartesian,
    linear_index

abstract type AbstractGridVector{G<:Grid} <: LinearCombination end

struct GridVector{T <: Grid} <: AbstractGridVector{T}
    grid::T
    _coefficients::SVector{3, Int}
    ket::Bool
end

make_grid_vector(g::T, _coefficients::AbstractVector) where {T<:Grid} =
    GridVector{T}(g, _coefficients, true)

grid(grid_vector::AbstractGridVector) = grid_vector.grid
set_grid(grid_vector::AbstractGridVector, new_grid) = @set grid_vector.grid = new_grid

miller_to_standard(grid_vector::AbstractGridVector, offsets::Tuple) =
    miller_to_standard(size(grid(grid_vector)), tuple(coefficients(grid_vector)...), offsets)

# miller_to_standard(grid_vector::AbstractGridVector{T}, offsets::AbstractGridVector{T}) where {T} =
#     miller_to_standard(grid_vector, coefficients(offsets))


"""
    wrapped(grid_vector)

The grid vector modulo the domain of the grid.

For grid vectors within the grid domain, this just gives the _coefficients.
For those outside the grid domain, this move the vector by multiples of 
the grid domain to bring it into the grid domain.
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
This gives the number of times we have to translate a vector by a grid domain to
bring them into the grid domain.
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

# function overflow(grid_vector::AbstractGridVector)
#     c = coefficients(grid_vector)
#     d = domain(grid(grid_vector))
#     result = zeros(Int, 3)
#     for i = 1:3
#         (l, r) = d[i]
#         result[i] = fld((c[i] - l), r - l + 1)
#     end
#     return result
# end


"""
The gives a image of the grid_vector within the grid domain.
"""
function reset_overflow(grid_vector::AbstractGridVector) 
    make_grid_vector(grid(grid_vector), wrapped(grid_vector))
end

"""
This moves the grid_vector by some multiples of the grid domain.
"""
function add_overflow(grid_vector::AbstractGridVector, overflow) 
    new_coefficients = map((o, c, s)->c + o * s, overflow, coefficients(grid_vector), SVector(size(grid(grid_vector))))
    make_grid_vector(grid(grid_vector), new_coefficients)
end

"""
Whether there is any overflow.
"""
has_overflow(grid_vector::AbstractGridVector) =
    any(f -> f != 0, overflow(grid_vector))


"""
Return the basis as they are internally stored (as kets).
"""
_basis(grid_vector::AbstractGridVector) = basis(grid(grid_vector))

function add(grid_vector_1::T, grid_vector_2::T) where T <: AbstractGridVector
    # ket(grid_vector_1) == ket(grid_vector_2) || error("Adding a bra to a ket.")
    grid(grid_vector_1) == grid(grid_vector_2) || error("Grid vectors defined on different grid.")
    new_coefficients = coefficients(grid_vector_1) + coefficients(grid_vector_2)
    return T(grid(grid_vector_1), new_coefficients, ket(grid_vector_1))
end

negate(grid_vector_1::T) where T <: AbstractGridVector =
    T(grid(grid_vector_1), -coefficients(grid_vector_1), ket(grid_vector_1))
mul(s::Int, l1::T) where T <: AbstractGridVector = T(grid(l1), s * coefficients(l1), ket(l1))
# minus(grid_vector_1::T, grid_vector_2::T) where T <: AbstractGridVector = add(grid_vector_1, negate(grid_vector_2))

function Base.show(io::IO, grid_vector::AbstractGridVector)
    print(io, "$(typeof(grid_vector)):\n")
    # print(io, "grid: $(typeof(grid(grid_vector)))\n" |> indent)
    print(io, indent("_coefficients: $(coefficients(grid_vector))") * "\n")
end

Base.:(==)(grid_vector_1::T, grid_vector_2::T) where T <: AbstractGridVector =
    coefficients(grid_vector_1) == coefficients(grid_vector_2) # && grid(grid_vector_1) == grid(grid_vector_2) 

Base.hash(grid_vector::AbstractGridVector)::UInt = linear_index(grid_vector) |> abs

"""
    cartesian(grid_vector)

The Cartesian coordinates of a grid vector.
"""
cartesian(grid_vector::AbstractGridVector)::Vector{Number} = let g = grid(grid_vector)
    basis_matrix(g) * (coefficients(grid_vector) + SVector(center(g)))
end
# cartesian(grid_vec::AbstractGridVector)::Vector{Number} = basis_transform(
#     coefficients(grid_vec), basis(grid_vec), CARTESIAN_BASIS)

LinearAlgebra.norm(grid_vector::AbstractGridVector) = norm(cartesian(grid_vector))

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
