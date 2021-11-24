export OnGrid,
    grid,
    grid!,
    translate,
    translate!,
    translation,
    element,
    element!,
    elements,
    elements!,
    standardize,
    wtp_normalize,
    wtp_normalize!

"""
A function defined on a grid. 
Examples include orbitals and bands.

Centering Convention:

An OnGrid complies to the centering convention if it is defined
on a grid centered at the origin. This means the domain is
from -N to N-1 for an even grid, and -N to N for an odd grid.
"""
abstract type OnGrid{G<:Grid} end

"""
The grid on which the OnGrid object is defined.
"""
grid(on_grid::OnGrid)::Grid = on_grid.grid
grid!(on_grid::OnGrid, new_grid) = on_grid.grid = new_grid

"""
Access to the elements of on_grid as they are internally stored.
"""
element(on_grid::OnGrid, indices...) = on_grid.elements[indices...]
element!(on_grid::OnGrid, val, indices...) = on_grid.elements[indices...] = val

elements(on_grid::OnGrid) = on_grid.elements
elements!(on_grid::OnGrid, new_elements) = on_grid.elements = new_elements



miller_to_standard(grid_vector::GridVector, offsets::AbstractVector{Int64}) =
    miller_to_standard(size(grid(grid_vector)), coefficients(grid_vector), offsets)

miller_to_standard(grid_vector::GridVector{T}, offsets::GridVector{T}) where {T} =
    miller_to_standard(grid_vector, coefficients(offsets))

"""
Indexing an OnGrid object with a grid vector gives the 
element on the corresponding grid point.
"""
function Base.getindex(on_grid::OnGrid, grid_vector::GridVector)
    !has_overflow(grid_vector) || error("overflow: $(grid_vector)\n on \n$(grid(on_grid))")
    grid(grid_vector) == grid(on_grid) || error("mismatching grid")
    offsets = translation(on_grid)
    indices = miller_to_standard(grid_vector, offsets)
    return element(on_grid, indices...)
end

function Base.getindex(on_grid::OnGrid, grid_vector_array::AbstractArray{<:GridVector})
    # TODO: Implement error checking.
    offsets = translation(on_grid)
    index_array = (v->miller_to_standard(v, offsets)).(grid_vector_array)
    return (I->element(on_grid, I...)).(index_array)
end

function Base.setindex!(on_grid::OnGrid, value, grid_vector::GridVector)
    !has_overflow(grid_vector) || error("overflow: $(grid_vector)\n on \n$(grid(on_grid))")
    grid(grid_vector) == grid(on_grid) || error("mismatching grid")
    offsets = translation(on_grid)
    indices = miller_to_standard(grid_vector, offsets)
    element!(on_grid, value, indices...)
end

function Base.setindex!(on_grid::OnGrid, value_array::AbstractArray, grid_vector_array::AbstractArray{<:GridVector})
    offsets = translation(on_grid)
    index_array = (v->miller_to_standard(v, offsets)).(grid_vector_array)
    map((v, i)->element!(on_grid, v, i...), value_array, index_array)
end

"""
Transalate/move on_grid by amount. 
This is done virtually by shifting the domain of the underlying grid.
"""
function translate(on_grid::OnGrid{T}, amount::GridVector{T}) where {T<:Grid}
    g = grid(on_grid)
    new_g = @set g.domain =
        Tuple((d[1] + t, d[2] + t) for (d, t) in zip(domain(g), coefficients(amount)))
    @set on_grid.grid = new_g
end

function translate!(on_grid::OnGrid{T}, amount::GridVector{T}) where {T<:Grid}
    g = grid(on_grid)
    new_g = @set g.domain =
        Tuple((d[1] + t, d[2] + t) for (d, t) in zip(domain(g), coefficients(amount)))
    grid!(on_grid, new_g)
end


function Base.string(on_grid::OnGrid)
    type_str = "type: $(typeof(on_grid))"
    grid_str = "grid:\n$(indent(repr(grid(on_grid))))"
    translation_str = "translation: $(translation(on_grid))"
    element_str = "element_type:\n$(indent(repr(Base.return_types(element, (typeof(on_grid), Integer, Integer, Integer))[1])))"
    return join([type_str, grid_str, translation_str, element_str], "\n")
end

"""
    translation(on_grid)

The translation of an OnGrid object is the center of the underlying grid as a
grid vector. 
"""
translation(on_grid::OnGrid) = let g = grid(on_grid)
    grid_vector_constructor(g, center(g))
end

"""
Normalize on_grid to unit norm. 
"""
function wtp_normalize(on_grid::OnGrid)
    elems = elements(on_grid)
    normalization_factor = norm(elems)
    @set on_grid.elements = elems / normalization_factor
end

"""
Same as normalize, but may modify the argument if mutable.
"""
function wtp_normalize!(on_grid::OnGrid)
    elems = elements(on_grid)
    normalization_factor = norm(elems)
    # println(normalization_factor)
    elements!(on_grid, elems / normalization_factor)
    return on_grid
end

"""

Standardize the representation of an OnGrid object.
The resulting object will be defined from -N+1 (N) to N.
Values outside the grid will be wrapped around.

                     2 ... N-1
-N+1 -N+2 ...  0  1  2 ... N-1 N

This operation does modify the underlying array.
This default implementation is slow.
"""
function standardize(on_grid::OnGrid)
    amount = translation(on_grid)
    non_zeros = filter(n -> n != 0, coefficients(amount))
    length(non_zeros) == 0 && return on_grid

    new_on_grid = translate(deepcopy(on_grid), -amount)
    for k in grid(new_on_grid)
        new_on_grid[k] = on_grid[k]
    end
    return new_on_grid
end

function Base.:(>>)(on_grid::OnGrid, translation::AbstractVector{<: Number})
    standardize(translate(on_grid, grid(on_grid)[translation...]))
end

struct SimpleFunctionOnGrid{T} <: OnGrid{T}
    grid::T
    elements::AbstractArray
    ket::Bool
end

function Base.map(f::Function, grid::T) where T <: Grid
    raw_elements = map(f, array(grid))
    shifted_elements = circshift(raw_elements, mins(grid))
    result = SimpleFunctionOnGrid{T}(grid, shifted_elements, true)
    return result
end

function resemble(on_grid, ::Type{SimpleFunctionOnGrid{T}}) where T
    g = grid(on_grid)
    SimpleFunctionOnGrid(g, zeros(T, size(g)), ket(g))
end

function add(o_1::T, o_2::S) where {T <: OnGrid, S <: OnGrid}
    grid(o_1) == grid(o_2) || error("Mismatching Grids.")
    ket(o_1) == ket(o_2) || error("Adding bra to ket.")
    o_3 = resemble(o_2, T)
    elements!(o_3, elements(o_1) + elements(o_2))
    return o_3
end

function negate(o_1::T) where T <: OnGrid
    o_2 = resemble(o_1, T)
    elements!(o_2, -elements(o_1))
    return o_2
end

function minus(o_1::OnGrid, o_2::OnGrid)
    add(o_1, negate(o_2))
end

function mul(o_1::T, o_2::S) where {T <: OnGrid, S <: OnGrid}
    grid(o_1) == grid(o_2) || error("Mismatching Grids.")
    ket(o_1) == ket(o_2) || error("elementwise product cannot take a bra and a ket.")
    o_3 = resemble(o_2, S)
    elements!(o_3, elements(o_1) .* elements(o_2))
    return o_3
end

function mul(scalar::Number, o_1::OnGrid)
    o_2 = resemble(o_1, T)
    elements!(o_2, scalar * elements(o_1))
    return o_2
end

Base.:+(o_1::OnGrid, o_2::OnGrid) = add(o_1, o_2)
Base.:-(o_1::OnGrid) = negate(o_1)
Base.:-(o_1::OnGrid, o_2::OnGrid) = minus(o_1, o_2)
Base.:*(scalar, o_1::OnGrid) = mul(scalar, o_1)
Base.:*(vector::AbstractVector, o_1::OnGrid) = [mul(s, o_1) for s in vector]
Base.:*(o_1::OnGrid, scalar) = mul(scalar, o_1)
Base.:*(o_1::OnGrid, vector::AbstractVector) = mul(vector, o_1)

Base.:*(o_1::OnGrid, o_2::OnGrid) = ket(o_1) == ket(o_2) ? mul(o_1, o_2) : braket(o_1, o_2)
Base.adjoint(o_1::OnGrid) = dagger(o_1)
