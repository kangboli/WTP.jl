export kpoint,
    kpoint!,
    i_kpoint,
    i_kpoint!,
    cache,
    cache!,
    index_band,
    compare,
    fft,
    AbstractUnkOrbital,
    KPoint,
    UnkBasisOrbital,
    UnkOrbital,
    orthonormal,
    ifft

"""
An UnkOrbital represents a u_{nk}(x) function. 

Such a function is also a basis vector and a linear combination.
This is where the single subtyping system of Julia gives me migraine.
"""
abstract type AbstractUnkOrbital{T} <: OnGrid{T} end
const KPoint = AbstractGridVector{<:BrillouinZone}

Base.adjoint(orbital::AbstractUnkOrbital) = dagger(orbital)
i_kpoint(orbital::AbstractUnkOrbital) = haskey(orbital.meta, :i_kpoint) ? orbital.meta[:i_kpoint] : nothing
i_kpoint!(orbital::AbstractUnkOrbital, i) = orbital.meta[:i_kpoint] = i
cache(orbital::AbstractUnkOrbital) = haskey(orbital.meta, :cache) ? orbital.meta[:cache] : nothing
cache!(orbital::AbstractUnkOrbital, c) = orbital.meta[:cache] = c

kpoint(orbital) = orbital.kpoint
function kpoint!(orbital, new_kpoint::KPoint) 
     orbital.kpoint = new_kpoint
     return orbital
end
index_band(orbital) = orbital.index_band

function braket(o_1::OnGrid, o_2::OnGrid)
    !ket(o_1) && ket(o_2) || error("braket requires a bra and a ket.")
    # cache lookup.
    # This cache turned out to be very slow.
    # c = cache(o_1)
    # if c !== nothing 
    #     result = c[kpoint(o_1), kpoint(o_2)]
    #     result !== nothing && return result[index_band(o_1), index_band(o_2)]
    #     println("no cache entry for $(kpoint(o_1)) \n and \n $(kpoint(o_2))")
    # end
    translation(o_1) == translation(o_2) || error("orbitals not aligned: $(translation(o_1))\n $(translation(o_2))")
    v_1 = reshape(elements(o_1), length(grid(o_1)))
    v_2 = reshape(elements(o_2), length(grid(o_2)))
    return transpose(v_1) * v_2
end

# Base.:*(o_1::Any, o_2::AbstractUnkOrbital) = braket(o_1, o_2)
# Base.:*(o_1::AbstractUnkOrbital, o_2::Any) = braket(o_1, o_2)

mutable struct UnkBasisOrbital{T} <: AbstractUnkOrbital{T}
    grid::T
    elements::AbstractArray
    kpoint::AbstractGridVector{<:BrillouinZone}
    index_band::Integer
    ket::Bool
    meta::Dict{Symbol,Any}
end

function UnkBasisOrbital(
    grid::T,
    elements::AbstractArray,
    kpoint::AbstractGridVector{<:BrillouinZone},
    index_band::Integer,
) where {T} 
    orbital = UnkBasisOrbital{T}(grid, elements, kpoint, index_band, true, Dict()) 
    return orbital
end 


## Override the functions of Basis.

function dagger(orbital::UnkBasisOrbital)
    orbital = @set orbital.ket = !orbital.ket
    orbital = @set orbital.elements = conj.(elements(orbital))
    return orbital
end

function dagger!(orbital::UnkBasisOrbital)
    orbital.elements = conj.(elements(orbital))
    orbital.ket = !orbital.ket
end

function resemble(orbital::UnkBasisOrbital{S}, ::Type{T}, new_elements=nothing) where {S <: Grid,  T <: Grid}
    g = grid(orbital)
    S == dual_grid(T) && (g = transform_grid(g))
    if new_elements === nothing 
       new_elements = zeros(eltype(elements(on_grid)), size(g))
    end
    UnkBasisOrbital(g, new_elements, kpoint(orbital), index_band(orbital)) 
end


# """
# Fast Fourier Transform of an orbital.

# This can be, and should be, parallelized with PencilFFT.jl
# """
# function fft(orbital::UnkBasisOrbital{T}) where T <:HomeCell
#     new_elements = FFTW.fft(elements(orbital))
#     new_grid = transform_grid(grid(orbital))

#     return UnkBasisOrbital{dual_grid(T)}(
#         new_grid,
#         new_elements,
#         kpoint(orbital),
#         index_band(orbital),
#         ket(orbital),
#         Dict(),
#     ) |> wtp_normalize!
# end

# function ifft(orbital::UnkBasisOrbital{T}) where T<:ReciprocalLattice
#     new_elements = FFTW.ifft(elements(orbital))
#     new_grid = transform_grid(grid(orbital))

#     return UnkBasisOrbital{dual_grid(T)}(
#         new_grid,
#         new_elements,
#         kpoint(orbital),
#         index_band(orbital),
#         ket(orbital),
#         Dict(),
#     ) |> wtp_normalize!
# end


function Base.show(io::IO, orbital::UnkBasisOrbital)
    ket(orbital) ? print(io, "ket\n") : print(io, "bra\n")
    print(io, "grid:\n$(indent(repr(grid(orbital))))" * "\n")
    print(io, "kpoint:\n$(indent(repr(kpoint(orbital))))\n")
    print(io, "band:\n$(indent(repr(index_band(orbital))))")
end

"""
A linear combination of UnkBasisOrbitals.

This also should be a subtype of OnGrid and Basis, but 
we are crippled by the single subtyping system.
"""
mutable struct UnkOrbital <: LinearCombination
    _coefficients::Vector{Number}
    basis::Any
    ket::Bool
    orthonormal::Bool
end

kpoint(orbital::UnkOrbital) =
    isempty(orbital.basis) ? nothing : kpoint(orbital.basis[1])
function kpoint!(orbital::UnkOrbital, new_kpoint::KPoint)
    (b -> kpoint!(b, new_kpoint)).(orbital.basis)
end
index_band(orbital::UnkOrbital) =
    isempty(orbital.basis) ? nothing : index_band(orbital.basis[1])
orthonormal(orbital::UnkOrbital) = orbital.orthonormal

UnkOrbital(orbital, orthonormal = true) = UnkOrbital(
    [1.0],
    [ket(orbital) ? orbital : dagger(orbital)],
    ket(orbital),
    orthonormal,
)

## Functions from Basis

function dagger(orbital::UnkOrbital)
    UnkOrbital(
        conj.(coefficients(orbital)),
        _basis(orbital),
        !ket(orbital),
        orthonormal(orbital),
    )
end
function dagger!(orbital::UnkOrbital)
    orbital._coefficients = conj.(coefficients(orbital))
    orbital.ket = !orbital.ket
end

function braket(o_1::UnkOrbital, o_2::UnkOrbital)
    !ket(o_1) && ket(o_2) || error("braket requires a bra and a ket.")
    o_1.basis === o_2.basis && orthonormal(o_1) && return transpose(coefficients(o_1)) * coefficients(o_2)

    return invoke(braket, Tuple{LinearCombination,LinearCombination}, o_1, o_2)
end

braket(o_1::UnkBasisOrbital, o_2::UnkOrbital) = braket(UnkOrbital(o_1), o_2)
braket(o_1::UnkOrbital, o_2::UnkBasisOrbital) = braket(o_1, UnkOrbital(o_2))


# Functions for LinearCombination

function LinearAlgebra.zeros(UnkOrbital, dims...)
    fill(UnkOrbital([], [], true, true), dims...)
end

function LinearAlgebra.zero(::UnkOrbital)
    UnkOrbital([], [], true, true)
end

function add(o_1::UnkOrbital, o_2::UnkOrbital, mutual_orthonormal = false)
    ket(o_1) == ket(o_2) || error("adding a bra to a ket.")
    o_1.basis === o_2.basis &&
        return UnkOrbital(coefficients(o_1) + coefficients(o_2), _basis(o_1), ket(o_1), orthonormal(o_1))

    coefficients_dict = OrderedDict()
    for (k, v) in zip([_basis(o_1)..., _basis(o_2)...], [coefficients(o_1)..., coefficients(o_2)...])
        coefficients_dict[k] = haskey(coefficients_dict, k) ? coefficients_dict[k] + v : v
    end

    # Basis orthogonality is not checked.

    UnkOrbital(
        collect(values(coefficients_dict)),
        collect(keys(coefficients_dict)),
        ket(o_1),
        orthonormal(o_1) && orthonormal(o_2) && mutual_orthonormal,
    )
end

negate(orbital::UnkOrbital) = @set orbital._coefficients = -orbital._coefficients

mul(s::Number, orbital::UnkOrbital) =
    @set orbital._coefficients = s * orbital._coefficients


function Base.show(io::IO, orbital::UnkOrbital)
    ket(orbital) ? print(io, "ket\n") : print(io, "bra\n")
    print(io, "coefficients:\n$(indent(repr(coefficients(orbital))))" * "\n")
    print(io, "n_basis:\n$(indent(string(length(basis(orbital)))))" * "\n")
    # print(io, "kpoint:\n$(indent(repr(kpoint(orbital))))\n")
    # print(io, "band:\n$(indent(repr(index_band(orbital))))")
end

function compare(o_1, o_2)
    c_1 = coefficients(kpoint(o_1))
    c_2 = coefficients(kpoint(o_2))
    c_1 == c_2 && return index_band(o_1) < index_band(o_2)
    c_1[1] == c_2[1] && c_1[2] == c_2[2] && return c_1[3] < c_2[3]
    c_1[1] == c_2[1] && return c_1[2] < c_2[2]
    return c_1[1] < c_2[1]
end

consolidate(orbital::UnkBasisOrbital) = (grid(orbital), elements(orbital))

function consolidate(orbital::UnkOrbital)
    elements = nothing
    grid = nothing
    for (c, b) in zip(coefficients(orbital), basis(orbital))
        new_grid, new_elements = consolidate(b)    
        grid === nothing && (grid = new_grid)
        grid == new_grid || error("Grids Mismatching.")
        elements === nothing && (elements = zeros(ComplexFxx, size(new_grid)))
        elements += c * new_elements
    end
    return grid, elements
end

function UnkBasisOrbital(orbital::UnkOrbital)
    grid, elements = consolidate(orbital)
    UnkBasisOrbital( grid, elements, kpoint(orbital), index_band(orbital))
end
