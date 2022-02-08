export OnGrid,
    grid,
    # grid!,
    translate,
    # translate!,
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
    square_normalize,
    square_normalize!,
    sparsify,
    # sparsify!,
    expand,
    center_spread,
    compute_r2

"""
The orbitals are defined somewhat abstractly as functions on grids. The grid
points are often times either real space points or wave numbers, but more
generally they are sets of parameters of the a family of functions.

Centering Convention:

An OnGrid complies to the centering convention if it is defined
on a grid centered at the origin. This means the domain is
from -N to N-1 for an even grid, and -N to N for an odd grid.
"""
abstract type OnGrid{G<:Grid} end

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

Example

```jldoctest on_grid
julia> homecell = make_grid(HomeCell3D, CARTESIAN_BASIS, size_to_domain((4, 4, 4)));

julia> lattice = transform_grid(homecell);

julia> ϕ = map(r->exp(1im * (lattice[1, 0, 0]' * r)), homecell);

julia> ϕ[homecell[:, 0, 0]]
4-element Vector{ComplexF64}:
                  -1.0 - 1.2246467991473532e-16im
 6.123233995736766e-17 - 1.0im
                   1.0 + 0.0im
 6.123233995736766e-17 + 1.0im
```

This syntax is very simple and expressive, but keep in mind that this involves
compiling an anonymous function for each map, which costs as much as the mapping
it self. 
"""
function Base.map(f::Function, grid::T, threading=true) where T<:Grid
    grid_vectors = array(grid)
    raw_elements = threading ? Folds.map(f, grid_vectors, ThreadedEx()) : map(f, grid_vectors)
    shifted_elements = circshift(raw_elements, mins(grid))
    result = SimpleFunctionOnGrid{T}(grid, shifted_elements, true)
    return result
end


"""
    grid(on_grid)

The grid on which the OnGrid object is defined.
"""
grid(on_grid::OnGrid)::Grid = on_grid.grid
# grid!(on_grid::OnGrid, new_grid) = on_grid.grid = new_grid
"""
    set_grid(on_grid, new_grid)

Create a copy of `on_grid` with a `new_grid`.
"""
set_grid(on_grid::OnGrid, new_grid) = @set on_grid.grid = new_grid

"""
Access to the elements of on_grid as they are internally stored.
"""
element(on_grid::OnGrid, indices...) = on_grid.elements[indices...]
element!(on_grid::OnGrid, val, indices...) = on_grid.elements[indices...] = val

elements(on_grid::OnGrid) = on_grid.elements
elements!(on_grid::OnGrid, new_elements) = on_grid.elements = new_elements
set_elements(on_grid::OnGrid, new_elements) = @set on_grid.elements = new_elements


"""
    getindex(on_grid, grid_vector)

Indexing an `OnGrid` object with a grid vector gives the element on the
corresponding grid point. Can also write `on_grid[grid_vector]`

Example:

```jldoctest on_grid
julia> ϕ[homecell[0, 0, 0]]
0.0
```
"""
function Base.getindex(on_grid::OnGrid, grid_vector::AbstractGridVector)
    overflow_detection && has_overflow(grid_vector) && error("overflow: $(grid_vector)\n on \n$(grid(on_grid))")
    grid(grid_vector) == grid(on_grid) || error("mismatching grid")
    indices = miller_to_standard(grid_vector, center(grid(on_grid)))
    return element(on_grid, indices...)
end


"""
    getindex(on_grid, grid_vector_array)

Indexing an `OnGrid` object with an array of grid vectors gives an array of elements
with the same dimension.

Example:

```jldoctest on_grid
julia> ϕ[homecell[:, 0, -1:1]]
4×3 Matrix{Float64}:
 0.0  0.0  0.0
 0.0  0.0  0.0
 0.0  0.0  0.0
 0.0  1.0  0.0
```
"""
function Base.getindex(on_grid::OnGrid, grid_vector_array::AbstractArray{<:AbstractGridVector})
    # TODO: Implement error checking.
    index_array = (v -> miller_to_standard(v, center(grid(on_grid)))).(grid_vector_array)
    return (I -> element(on_grid, I...)).(index_array)
end

"""
    setindex!(on_grid, value, grid_vector)

Set an element of an `OnGrid` object identified by a grid vector.
Can also write `on_grid[grid_vector] = value`.

Example:

```jldoctest on_grid
julia> ϕ[homecell[0, 0, 0]] = 24;
julia> ϕ[homecell[0, 0, 0]]
24
julia> ϕ[homecell[0, 0, 0]] = 0; # remember to set it back for later examples.
```
"""
function Base.setindex!(on_grid::OnGrid, value, grid_vector::AbstractGridVector)
    overflow_detection && has_overflow(grid_vector) && error("overflow: $(grid_vector)\n on \n$(grid(on_grid))")
    grid(grid_vector) == grid(on_grid) || error("mismatching grid")
    indices = miller_to_standard(grid_vector, center(grid(on_grid)))
    element!(on_grid, value, indices...)
end

"""
    setindex!(on_grid, value, grid_vector)

Set an array of elements of an `OnGrid` object identified by an array of grid
vectors with the same dimension.
Can also write `on_grid[grid_vector] = value`.

Example:

```jldoctest on_grid
julia> ϕ[homecell[:, 0, 0]] = [4, 3, 2, 1];
julia> ϕ[homecell[:, 0, 0]]
24
julia> ϕ[homecell[:, 0, 0]] = [0, 0, 0, 1]; # remember to set it back for later examples.
```
"""
function Base.setindex!(on_grid::OnGrid, value_array::AbstractArray, grid_vector_array::AbstractArray{<:AbstractGridVector})
    index_array = (v -> miller_to_standard(v, center(grid(on_grid)))).(grid_vector_array)
    map((v, i) -> element!(on_grid, v, i...), value_array, index_array)
end

"""
Transalate/move on_grid by amount. 
This is done virtually by shifting the domain of the underlying grid.
"""
function translate(on_grid::OnGrid{T}, amount::AbstractVector) where {T<:Grid}
    g = grid(on_grid)
    new_g = set_domain(g, Tuple((d[1] + t, d[2] + t) for (d, t) in zip(domain(g), amount)))
    set_grid(on_grid, new_g)
end

# function translate!(on_grid::OnGrid{T}, amount::AbstractGridVector{T}) where {T<:Grid}
#     g = grid(on_grid)
#     new_g = @set g.domain =
#         Tuple((d[1] + t, d[2] + t) for (d, t) in zip(domain(g), coefficients(amount)))
#     grid!(on_grid, new_g)
# end


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
translation(on_grid::OnGrid) = center(grid(on_grid))

"""
    norm(on_grid)

The square norm of a function on a grid.

Example:

```jldoctest on_grid
julia> norm(ϕ)
8.0
```
"""
function LinearAlgebra.norm(on_grid::OnGrid)
    return on_grid |> elements |> norm
end

"""
    square_normalize(on_grid)

Normalize on_grid to unit square norm. 

Example:

```jldoctest on_grid
julia> ϕ_2 = square_normalize(ϕ);
julia> norm(ϕ_2)
1.0
```
"""
function square_normalize(on_grid::OnGrid)
    elems = elements(on_grid)
    normalization_factor = norm(elems)
    @set on_grid.elements = elems / normalization_factor
end

"""
    square_normalize!(on_grid)

Same as `square_normalize`, but in place.

Example:

```jldoctest on_grid
julia> square_normalize!(ϕ);
julia> norm(ϕ)
1.0
```
"""
function square_normalize!(on_grid::OnGrid)
    normalization_factor = norm(elements(on_grid))
    rdiv!(elements(on_grid), normalization_factor)
    return on_grid
end

"""
    fft(on_grid)

Fast Fourier Transform of an `OnGrid` object.

Example:

```jldoctest on_grid
julia> ϕ̃ = fft(ϕ);
julia> ϕ̃[lattice[:, 0, 0]]
4-element Vector{ComplexF64}:
 -3.061616997868383e-17 - 3.061616997868383e-17im
                    0.0 + 3.061616997868383e-17im
  3.061616997868383e-17 - 3.061616997868383e-17im
                    1.0 + 3.061616997868383e-17im

julia> typeof(ψ)
WTP.SimpleFunctionOnGrid{ReciprocalLattice3D}
```
"""
function fft(on_grid::OnGrid{T}, normalize = true) where {T<:HomeCell}
    new_elements = FFTW.fft(elements(on_grid))
    unnormalized = resemble(on_grid, dual_grid(T), new_elements)
    return normalize ? unnormalized |> square_normalize! : unnormalized
end

"""
    _fft!(on_grid)

Generally you shouldn't use this method since it leaves `on_grid`
in an invalid state. Use `@fft!` instead to set `on_grid` to `nothing`.

In place Fast Fourier Transform (FFT) of an `OnGrid` object.
The memory of the argument will be repurposed for its FFT.
"""
function _fft!(on_grid::OnGrid{T}, normalize = true) where {T<:HomeCell}
    FFTW.fft!(elements(on_grid))
    unnormalized = resemble(on_grid, dual_grid(T), elements(on_grid))
    return normalize ? unnormalized |> square_normalize! : unnormalized
end

"""
    @fft!(real_orbital)

Perform a in-place Fourier Transform. 
`real_orbital` will be set to `nothing`.
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

Example:

```jldoctest on_grid
julia> ϕ = ifft(ϕ̃)

julia> ϕ[homecell[:, 0, 0]]
4-element Vector{ComplexF64}:
                -0.125 - 1.5308084989341915e-17im
 7.654042494670958e-18 - 0.125im
                 0.125 + 0.0im
 7.654042494670958e-18 + 0.125im
julia> typeof(ϕ)
WTP.SimpleFunctionOnGrid{HomeCell3D}
```
"""
function ifft(on_grid::OnGrid{T}, normalize = true) where {T<:ReciprocalLattice}
    new_elements = FFTW.ifft(elements(on_grid))
    unnormalized = resemble(on_grid, dual_grid(T), new_elements)
    return normalize ? unnormalized |> square_normalize! : unnormalized
end

"""
    _ifft!(orbital)

Generally you shouldn't use this method since it leaves `on_grid`
in an invalid state. Use `@ifft!` instead to set `on_grid` to `nothing`.

In place Inverse Fast Fourier Transform of an `OnGrid` object.
"""
function _ifft!(on_grid::OnGrid{T}, normalize = true) where {T<:ReciprocalLattice}
    FFTW.ifft!(elements(on_grid))
    unnormalized = resemble(on_grid, dual_grid(T), elements(on_grid))
    return normalize ? unnormalized |> square_normalize! : unnormalized
end

"""
    @ifft!(reciprocal_orbital)

Perform a in-place inverse Fourier Transform. 
`reciprocal_orbital` will be set to `nothing`.
"""
macro ifft!(reciprocal_orbital)
    esc(Expr(:block, Expr(
            :(=), :real_orbital,
            Expr(:call, Expr(Symbol("."), :WTP, QuoteNode(:_ifft!,)), reciprocal_orbital)
        ),
        Expr(:(=), reciprocal_orbital, :nothing),
        :real_orbital))
end

"""
    standardize(orbital)

Standardize the representation of an OnGrid object.
The resulting object will be defined from -N+1 (N) to N.
Values outside the grid will be wrapped around.

                     2 ... N-1
-N+1 -N+2 ...  0  1  2 ... N-1 N
"""
function standardize(orbital::OnGrid)
    amount = translation(orbital)
    non_zeros = filter(n -> n != 0, amount)
    length(non_zeros) == 0 && return orbital

    new_orbital = translate(orbital, -[amount...])
    elements!(new_orbital, zeros(ComplexFxx, size(grid(orbital))))
    circshift!(elements(new_orbital), elements(orbital), amount)
    return new_orbital
end

function Base.:(>>)(on_grid::OnGrid, translation::AbstractVector{<:Number})
    standardize(translate(on_grid, translation))
end

function resemble(on_grid::SimpleFunctionOnGrid{S}, ::Type{T}, new_elements = nothing) where {S<:Grid,T<:Grid}
    g = grid(on_grid)
    S == dual_grid(T) && (g = transform_grid(g))
    if new_elements === nothing
        new_elements = zeros(eltype(elements(on_grid)), size(g))
    end
    SimpleFunctionOnGrid(g, new_elements, ket(on_grid))
end

"""
    add(o_1, o_2)
    
Can also write `o_1 + o_2`
"""
function add(o_1::OnGrid{T}, o_2::OnGrid{T}) where {T<:Grid}
    grid(o_1) == grid(o_2) || error("Mismatching Grids.")
    # ket(o_1) == ket(o_2) || error("Adding bra to ket.")
    o_3 = resemble(o_2, T, elements(o_1) + elements(o_2))
    # elements!(o_3, )
    # o_3 = set_elements(o_2, elements(o_1) + elements(o_2))
    return o_3
end

"""
    negate(o_1)

Can also write `-o_1`.
"""
function negate(o_1::OnGrid{T}) where {T<:Grid}
    o_2 = resemble(o_1, T, -elements(o_1))
    # elements!(o_2, )
    return o_2
end


"""
    minus(o_1, o_2)

Can also write `o_1 - o_2`.
"""
function minus(o_1::OnGrid{T}, o_2::OnGrid{T}) where T <: Grid
    add(o_1, negate(o_2))
end


"""
    mul(o_1, o_2)

Elementwise product. Can also write `o_1 * o_2` when both of them are `bra` or `ket`.
"""
function mul(o_1::OnGrid{T}, o_2::OnGrid{S}) where {T<:Grid,S<:Grid}
    grid(o_1) == grid(o_2) || error("Mismatching Grids.")
    ket(o_1) == ket(o_2) || error("elementwise product cannot take a bra and a ket.")
    # o_3 = resemble(o_2, S)
    # elements!(o_3, elements(o_1) .* elements(o_2))
    o_3 = set_elements(o_2, elements(o_1) .* elements(o_2))
    return o_3
end

"""
    mul(scalar, o_1)

Scalar multiply. Can also write `scalar * o_1` or `o_1 * scalar`.
"""
function mul(scalar::Number, o_1::OnGrid{T}) where {T}
    o_2 = resemble(o_1, T, scalar * elements(o_1))
    # elements!(o_2, )
    return o_2
end

"""
    abs2(o_1)

Elementwise `abs2`. Handy for computing the density.

Example:

```jldoctest on_grid
julia> sum(elements(abs2(ϕ)))
1.0
```
"""
function Base.abs2(o_1::OnGrid{T}) where {T}
    o_2 = resemble(o_1, T, abs2.(elements(o_1)))
    return o_2
end

"""
    braket(o_1, o_2)

⟨o_1 | o_2⟩. It does not matter whether `o_1` or `o_2` is a ket. `braket` will
do the right thing. Can also write `o_1' * o_2`, where `o_1'` has to be a bra
and `o_2` must be a ket.

Example:

```jldoctest on_grid
julia> braket(ϕ, ϕ)
1.0 + 0.0im
julia> ψ' * ψ
1.0 + 6.123233995736766e-17im
```
"""
function braket(o_1::OnGrid, o_2::OnGrid)
    # !ket(o_1) && ket(o_2) || error("braket requires a bra and a ket.")
    # translation(o_1) == translation(o_2) || error("orbitals not aligned: $(translation(o_1))\n $(translation(o_2))")
    v_1 = reshape(elements(o_1), length(grid(o_1)))
    v_2 = reshape(elements(o_2), length(grid(o_2)))
    return (!ket(o_1) ? transpose(v_1) : adjoint(v_1)) * (ket(o_2) ? v_2 : conj(v_2))
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


function html(on_grid::OnGrid)
    return "<ul>
    <li>Type: $(typeof(on_grid))</li>
    <li>Grid: $(html(grid(on_grid)))</li>
    </ul>
    "
end

Base.show(io::IO, ::MIME"text/html", on_grid::OnGrid) = println(io, html(on_grid))

"""
    r2(homecell)

Compute ``r^2`` on a homecell together with its Fourier transform.

Example:

```jldoctest on_grid
julia> r2, r̃2 = compute_r2(homecell);
julia> r2[homecell[2, 0, 0]]
4.0
julia> r̃2[grid(r̃2)[0, 0, 0]]
288.0 + 0.0im
```
"""
function compute_r2(g::HomeCell)
    r_coordinates = cartesian.(g(1:length(g)))
    r2 = norm.(r_coordinates).^2
    r2 = SimpleFunctionOnGrid(g, reshape(r2, size(g)...), true)
    return r2, fft(r2, false)
end

"""
    center_spread(õ, r̃2)

Compute the center and the spread. Here, `õ` should generally be 
fft of the density (instead of the orbital).

The algorithm can fail when the center is ill-defined

Example:

```jldoctest on_grid
julia> center_spread(fft(abs2(ϕ), false), r̃2)
([NaN, NaN, NaN], NaN)
julia> center_spread(fft(abs2(square_normalize(map(r->r==homecell[0, 0, 0], homecell)))), r̃2)
([0.0, -0.0, -0.0], 8.326672684688674e-17)
```
"""
function center_spread(õ::OnGrid{T}, r̃2::OnGrid{T}) where {T<:ReciprocalLattice}
    convolved = ifft(õ * r̃2, false)
    elements!(convolved, abs.(elements(convolved)))
    linear_index_min = argmin(reshape(elements(convolved), length(grid(convolved))))
    r_min = grid(convolved)(linear_index_min)
    return quadratic_fit(convolved, r_min)
end



"""
A seven-points fitting method for the convolution. 

f(r) = (r-b)ᵀ D (r - b) + c

D is diagonal because r² = x² + y² + z² is additively separable.
This is:

 c + ∑ᵢ₌₁³ bᵢ² dᵢ - 2 bᵢ dᵢ rᵢ + dᵢ rᵢ² = y

[1, rᵢ, rᵢ²] [c + ∑ bᵢ² dᵢ, -2 bᵢ dᵢ, dᵢ]' = y
"""
function quadratic_fit(o::OnGrid{T}, fitting_center::GridVector{T}) where {T<:Grid}
    g = grid(o)
    fitting_points = reset_overflow.((v -> fitting_center + v).([
        g[1, 0, 0], g[-1, 0, 0], g[0, 1, 0],
        g[0, -1, 0], g[0, 0, 1], g[0, 0, -1], g[0, 0, 0]]))

    A = vcat(map(r -> [1, cartesian(r)..., (cartesian(r) .^ 2)...]', fitting_points)...)
    rhs = o[fitting_points]

    solution = A \ rhs
    d = solution[5:7]
    b = solution[2:4] ./ (-2d)
    minima = [1, b..., (b .^ 2)...]' * solution
    return b, minima
end

"""
A three-points fitting method for conjugate gradient.

f(λ) ≈ a (λ - b)² + c around λ = 1
f(λ) ≈ a λ² - 2ab λ + b² + c

| 0  0  1 | | a      | = | f(0) |
| 1  1  1 | | -2ab   |   | f(1) |
| 4  2  1 | | b² + c |   | f(2) |

```julia
quadratic_fit_1d(x->2(x-3)^2 + 3)

# output

(3.0, 3.0)
```
"""
function quadratic_fit_1d(f::Function)
    a, negative_2ab, _ = [.5 -1 .5; -1.5 2 -.5; 1.0 .0 .0] * f.(0:2)
    b = -negative_2ab / (2a)
    return b, f(b)
end


"""
    expand(on_grid, factors)

Expand (copy) an `OnGrid` object by `factors` along the respective directions.

Example:

```jldoctest on_grid
julia> ϕ̂ = expand(ϕ);

julia> grid(ϕ̂)
type: ReciprocalLattice3D
domain: ((-4, 3), (-4, 3), (-4, 3))
basis:
    ket: 1.571, 0.000, 0.000
    ket: 0.000, 1.571, 0.000
    ket: 0.000, 0.000, 1.571

julia> grid(ϕ)
type: ReciprocalLattice3D
domain: ((-2, 1), (-2, 1), (-2, 1))
basis:
    ket: 1.571, 0.000, 0.000
    ket: 0.000, 1.571, 0.000
    ket: 0.000, 0.000, 1.571
```
"""
function expand(on_grid::OnGrid, factors = [2, 2, 2])
    new_elements = repeat(elements(on_grid), factors...)
    g = grid(on_grid)
    # shift_amount = map((f, s)->isodd(f) ? 0 : -(s÷2), factors, size(g))
    # new_elements = circshift(new_elements, shift_amount)
    new_on_grid = set_elements(on_grid, new_elements)
    new_grid = expand(g, factors)
    return set_grid(new_on_grid, new_grid)
end

"""
    vectorize(o)

Reshape the elements into a column vector.
"""
vectorize(o::OnGrid) = reshape(elements(o), prod(size(grid(o))))

"""
    sparsify(on_grid, threshold=1e-16)

Store the elements of `on_grid` as a sparse array. 

PS: Quantum Espresso seems to have reinvented sparse arrays in trying to
store their orbitals (with a spherical energy cutoff) efficeintly.

Example: 

```jldoctest on_grid
julia> ϕ̃_2 = sparsify(ϕ̃, threshold=1e-16);
julia> vectorize(ϕ̃_2)
64-element SparseVector{Number, Int64} with 1 stored entry:
  [2 ]  =  1.0+3.06162e-17im
```
"""
function sparsify(on_grid::OnGrid; threshold=1e-16)
    # none_zero_indices = findall(!iszero, vectorize(on_grid))
    rounded = (n->abs(n) < threshold ? 0 : n).(vectorize(on_grid))
    sparse_vector = sparse(rounded)
    @set on_grid.elements = reshape(sparse_vector, size(elements(on_grid))...)
end

# function sparsify!(on_grid::OnGrid)
#     none_zero_indices = findall(!iszero, vectorize(on_grid))
#     elements!(on_grid, reshape(sparsevec(elements(on_grid)), size(elements(on_grid))...))
# end


"""
    on_grid >> grid_vector

Translate `on_grid` by `grid_vector`.

Example: 

```jldoctest on_grid
julia> ϕ̃ = ϕ̃ >> lattice[1, 0, 0];
julia> ϕ̃[lattice[:, 0, 0]] 
4-element Vector{ComplexF64}:
                    1.0 + 3.061616997868383e-17im
 -3.061616997868383e-17 - 3.061616997868383e-17im
                    0.0 + 3.061616997868383e-17im
  3.061616997868383e-17 - 3.061616997868383e-17im
julia> ϕ̃[lattice[:, 0, 0]]
4-element Vector{ComplexF64}:
-3.061616997868383e-17 - 3.061616997868383e-17im
                    0.0 + 3.061616997868383e-17im
3.061616997868383e-17 - 3.061616997868383e-17im
                    1.0 + 3.061616997868383e-17im
```
"""
function Base.:(>>)(on_grid::OnGrid{S}, grid_vector::AbstractGridVector{S}) where {S}
    on_grid >> coefficients(grid_vector)
end