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
    set_elements,
    fft,
    @fft!,
    _fft!,
    ifft,
    @ifft!,
    _ifft!,
    standardize,
    wtp_normalize,
    wtp_normalize!,
    wtp_sparse,
    wtp_sparse!,
    expand

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
set_grid(on_grid::OnGrid, new_grid) = @set on_grid.grid = new_grid

"""
Access to the elements of on_grid as they are internally stored.
"""
element(on_grid::OnGrid, indices...) = on_grid.elements[indices...]
element!(on_grid::OnGrid, val, indices...) = on_grid.elements[indices...] = val

elements(on_grid::OnGrid) = on_grid.elements
elements!(on_grid::OnGrid, new_elements) = on_grid.elements = new_elements
set_elements(on_grid::OnGrid, new_elements) = @set on_grid.elements = new_elements



miller_to_standard(grid_vector::AbstractGridVector, offsets::AbstractVector{Int64}) =
    miller_to_standard(size(grid(grid_vector)), coefficients(grid_vector), offsets)

miller_to_standard(grid_vector::AbstractGridVector{T}, offsets::AbstractGridVector{T}) where {T} =
    miller_to_standard(grid_vector, coefficients(offsets))

"""
Indexing an OnGrid object with a grid vector gives the 
element on the corresponding grid point.
"""
function Base.getindex(on_grid::OnGrid, grid_vector::AbstractGridVector)
    !has_overflow(grid_vector) || error("overflow: $(grid_vector)\n on \n$(grid(on_grid))")
    grid(grid_vector) == grid(on_grid) || error("mismatching grid")
    offsets = translation(on_grid)
    indices = miller_to_standard(grid_vector, offsets)
    return element(on_grid, indices...)
end

function Base.getindex(on_grid::OnGrid, grid_vector_array::AbstractArray{<:AbstractGridVector})
    # TODO: Implement error checking.
    offsets = translation(on_grid)
    index_array = (v->miller_to_standard(v, offsets)).(grid_vector_array)
    return (I->element(on_grid, I...)).(index_array)
end

function Base.setindex!(on_grid::OnGrid, value, grid_vector::AbstractGridVector)
    !has_overflow(grid_vector) || error("overflow: $(grid_vector)\n on \n$(grid(on_grid))")
    grid(grid_vector) == grid(on_grid) || error("mismatching grid")
    offsets = translation(on_grid)
    indices = miller_to_standard(grid_vector, offsets)
    element!(on_grid, value, indices...)
end

function Base.setindex!(on_grid::OnGrid, value_array::AbstractArray, grid_vector_array::AbstractArray{<:AbstractGridVector})
    offsets = translation(on_grid)
    index_array = (v->miller_to_standard(v, offsets)).(grid_vector_array)
    map((v, i)->element!(on_grid, v, i...), value_array, index_array)
end

"""
Transalate/move on_grid by amount. 
This is done virtually by shifting the domain of the underlying grid.
"""
function translate(on_grid::OnGrid{T}, amount::AbstractGridVector{T}) where {T<:Grid}
    g = grid(on_grid)
    new_g = @set g.domain =
        Tuple((d[1] + t, d[2] + t) for (d, t) in zip(domain(g), coefficients(amount)))
    @set on_grid.grid = new_g
end

function translate!(on_grid::OnGrid{T}, amount::AbstractGridVector{T}) where {T<:Grid}
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
    wtp_normalize(on_grid)

Normalize on_grid to unit norm. 
"""
function wtp_normalize(on_grid::OnGrid)
    elems = elements(on_grid)
    normalization_factor = norm(elems)
    @set on_grid.elements = elems / normalization_factor
end

"""
    wtp_normalize!(on_grid)

Same as normalize, but in place.
"""
function wtp_normalize!(on_grid::OnGrid)
    normalization_factor = norm(elements(on_grid))
    rdiv!(elements(on_grid), normalization_factor)
    return on_grid
end

function wtp_sparse(on_grid::OnGrid)
    @set on_grid.elements = reshape(sparsevec(elements(on_grid)), size(elements(on_grid))... )
end
function wtp_sparse!(on_grid::OnGrid)
    elements!(on_grid, reshape(sparsevec(elements(on_grid)), size(elements(on_grid))... ))
end

"""
Standardize the representation of an OnGrid object.
The resulting object will be defined from -N+1 (N) to N.
Values outside the grid will be wrapped around.

                     2 ... N-1
-N+1 -N+2 ...  0  1  2 ... N-1 N
"""
function standardize(orbital::OnGrid)
    amount = translation(orbital)
    non_zeros = filter(n -> n != 0, coefficients(amount))
    length(non_zeros) == 0 && return orbital

    new_orbital = translate(orbital, -amount)
    elements!(new_orbital, zeros(ComplexFxx, size(grid(orbital))))
    circshift!(elements(new_orbital), elements(orbital), Tuple(coefficients(amount)))
    return new_orbital
end

function Base.:(>>)(on_grid::OnGrid, translation::AbstractVector{<: Number})
    standardize(translate(on_grid, grid(on_grid)[translation...]))
end

function Base.:(>>)(on_grid::OnGrid{S}, grid_vector::AbstractGridVector{S}) where S
    on_grid >> coefficients(grid_vector)
end

mutable struct SimpleFunctionOnGrid{T} <: OnGrid{T}
    grid::T
    elements::AbstractArray
    ket::Bool
end

"""
    map(f, grid)

Map a *pure* function `f` onto every grid vector on `grid`.
The map is threaded, so an unpure `f` would probabbly be unsafe.

More likely than not, you would use the `do` syntax

```julia
map(homecell) do r
    exp(1im * k' * r)
end
```

This syntax is very simple and expressive, but keep in mind that this involves
compiling an anonymous function for each map, which costs as much as the mapping
it self. 
"""
function Base.map(f::Function, grid::T) where T <: Grid
    grid_vectors = array(grid)
    raw_elements = Folds.map(f, grid_vectors, ThreadedEx())
    shifted_elements = circshift(raw_elements, mins(grid))
    result = SimpleFunctionOnGrid{T}(grid, shifted_elements, true)
    return result
end

function resemble(on_grid::SimpleFunctionOnGrid{S}, ::Type{T}, new_elements=nothing) where {S <: Grid, T <:Grid}
    g = grid(on_grid)
    S == dual_grid(T) && (g = transform_grid(g))
    if new_elements === nothing 
       new_elements = zeros(eltype(elements(on_grid)), size(g))
    end
    SimpleFunctionOnGrid(g, new_elements, ket(on_grid))
end

function add(o_1::OnGrid{T}, o_2::OnGrid{T}) where {T <: Grid}
    grid(o_1) == grid(o_2) || error("Mismatching Grids.")
    # ket(o_1) == ket(o_2) || error("Adding bra to ket.")
    o_3 = resemble(o_2, T)
    elements!(o_3, elements(o_1) + elements(o_2))
    # o_3 = set_elements(o_2, elements(o_1) + elements(o_2))
    return o_3
end

function negate(o_1::OnGrid{T}) where T <: Grid
    o_2 = resemble(o_1, T)
    elements!(o_2, -elements(o_1))
    return o_2
end

function minus(o_1::OnGrid, o_2::OnGrid)
    add(o_1, negate(o_2))
end

function mul(o_1::OnGrid{T}, o_2::OnGrid{S}) where {T <: Grid, S <: Grid}
    grid(o_1) == grid(o_2) || error("Mismatching Grids.")
    ket(o_1) == ket(o_2) || error("elementwise product cannot take a bra and a ket.")
    # o_3 = resemble(o_2, S)
    # elements!(o_3, elements(o_1) .* elements(o_2))
    o_3 = set_elements(o_2, elements(o_1) .* elements(o_2))
    return o_3
end

function mul(scalar::Number, o_1::OnGrid{T}) where T
    o_2 = resemble(o_1, T)
    elements!(o_2, scalar * elements(o_1))
    return o_2
end

function Base.abs2(o_1::OnGrid{T}) where T
    o_2 = resemble(o_1, T, abs2.(elements(o_1)))
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


"""
    fft(on_grid)

Fast Fourier Transform of an `OnGrid` object.
"""
function fft(on_grid::OnGrid{T}) where T <: HomeCell
    new_elements = FFTW.fft(elements(on_grid))
    return resemble(on_grid, dual_grid(T), new_elements) |> wtp_normalize!
end

"""
    _fft!(on_grid)

Generally you shouldn't use this method since it leaves `on_grid`
in an invalid state. Use `@fft!` instead to set `on_grid` to `nothing`.

In place Fast Fourier Transform (FFT) of an `OnGrid` object.
The memory of the argument will be repurposed for its FFT.
"""
function _fft!(on_grid::OnGrid{T}) where T <: HomeCell 
    FFTW.fft!(elements(on_grid))
    return resemble(on_grid, dual_grid(T), elements(on_grid)) |> wtp_normalize!
end

"""
    @fft!(real_orbital)

Perform a in-place Fourier Transform. 
`real_orbital` will be set to `nothing`
"""
macro fft!(real_orbital)
    esc(Expr(:block, Expr(
        :(=), :reciprocal_orbital,
        Expr(:call, Expr(Symbol("."), :WTP, QuoteNode(:_fft!,)), real_orbital)
    ),
    Expr(:(=), real_orbital, :nothing),
    :reciprocal_orbital))
end

"""
    ifft(orbital)

Inverse Fast Fourier Transform of an `OnGrid` object.
The memory of the argument will be repurposed for its FFT.
"""
function ifft(on_grid::OnGrid{T}) where T <: ReciprocalLattice
    new_elements = FFTW.ifft(elements(on_grid))
    return resemble(on_grid, dual_grid(T), new_elements) |> wtp_normalize!
end

"""
    _ifft!(orbital)

Generally you shouldn't use this method since it leaves `on_grid`
in an invalid state. Use `@ifft!` instead to set `on_grid` to `nothing`.

In place Inverse Fast Fourier Transform of an `OnGrid` object.
"""
function _ifft!(on_grid::OnGrid{T}) where T <: ReciprocalLattice
    FFTW.ifft!(elements(on_grid))
    on_grid = resemble(on_grid, dual_grid(T), elements(on_grid)) |> wtp_normalize!
end

"""
    @ifft!(reciprocal_orbital)

Perform a in-place inverse Fourier Transform. 
`reciprocal_orbital` will be set to `nothing`
"""
macro ifft!(reciprocal_orbital) esc(Expr(:block, Expr(
        :(=), :real_orbital,
        Expr(:call, Expr(Symbol("."), :WTP, QuoteNode(:_ifft!,)), reciprocal_orbital)
    ),
    Expr(:(=), reciprocal_orbital, :nothing),
    :real_orbital))
end

function html(on_grid::OnGrid)
    return "<ul>
    <li>Type: $(typeof(on_grid))</li>
    <li>Grid: $(html(grid(on_grid)))</li>
    </ul>
    "
end

Base.show(io::IO, ::MIME"text/html", on_grid::OnGrid) = println(io, html(on_grid))

"""
    expand(on_grid, factors)

Expand (copy) an `OnGrid` object by `factors` along the respective directions.
"""
function expand(on_grid::OnGrid, factors=[2, 2, 2])
    new_elements = repeat(elements(on_grid), factors...)
    new_on_grid = set_elements(on_grid, new_elements)
    new_grid = expand(grid(on_grid), factors)
    return set_grid(new_on_grid, new_grid)
end

function vectorize(o::OnGrid)
    # wfc(o)!==nothing ? wfc(o).evc[:, index_band(o)] : 
    reshape(elements(o), prod(size(grid(o))))
end