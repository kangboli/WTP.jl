using SparseArrays

export NumericalOperator, integrate, indices, expand_storage

"""
The users should treat the integrals as an efficient lookup table.
They shouldn't have to worry about the details on how they are computed 
and stored. 

Typically, users have to compute the integrals and store them as a matrix. There
are two problems with this.

1. The mental models is integer indices and matrices instead of orbitals and
operators as they are Mathematically formulated.
2. Thinking in terms of matrices is complicated since the integrals are rarely simply 
matrices. It is often time in a sparse or compressed format.

We will discourage the users from dealing with raw matrices unless they must.
Integrals should be written as `ψ' | V | ϕ`. Accessing a "column" of `V` should
be `ϕ' | V | :`. The integrals should be computed, cached, and compressed
without user intervention until the stage of performance tweak. 
"""

mutable struct NumericalOperator{T} <: OnGrid{T}
    grid::T
    elements::AbstractArray{<:Number, <:Any}
    current_index::Integer
    storage::AbstractMatrix{ComplexFxx}
    valid::AbstractMatrix{Bool}
    indices::Dict{UInt64, <:Integer}
end

function NumericalOperator(grid::T, elements::AbstractArray{<:Number, <:Any}) where T <: Grid
    storage = zeros(ComplexFxx, (1, 1))
    valid = zeros(Bool, (1, 1))
    indices = Dict{UInt64, Integer}()
    NumericalOperator{T}(grid, elements, 1, storage, valid, indices)
end

function integrate(A::NumericalOperator{T}, o_1::OnGrid{T}, o_2::OnGrid{T}) where T <: Grid
    (!ket(o_1) && ket(o_2)) || error("A bra and a ket are required for integrals.")
    return sum(elements(o_1) .* elements(A) .* elements(o_2))
end


# Ugly hack with type Union due to the lack of multiple subtyping.
const Operator = Union{NumericalOperator}

cache(A::Operator) = A.cache
indices(A::Operator) = A.indices
current_index(A::Operator) = A.current_index
current_index!(A::Operator, new_index) = A.current_index = new_index
storage(A::Operator) = A.storage
storage!(A::Operator, new_storage) = A.storage = new_storage
valid(A::Operator) = A.valid
valid!(A::Operator, new_valid) = A.valid = new_valid

function expand_storage(A)
    current_size = size(storage(A), 1)
    new_size = current_size + max(current_size ÷ 2, 1)
    new_storage = zeros(ComplexFxx, (new_size, new_size))
    new_valid = zeros(Bool, (new_size, new_size))
    new_storage[1:current_size, 1:current_size] = storage(A)[:, :]
    new_valid[1:current_size, 1:current_size] = valid(A)[:, :]
    storage!(A, new_storage)
    valid!(A, new_valid)
end

function Base.getindex(A::Operator, o_1::OnGrid, o_2::OnGrid)
    if haskey(indices(A), objectid(o_1)) && haskey(indices(A), objectid(o_2)) 
        m, n = indices(A)[objectid(o_1)], indices(A)[objectid(o_2)]
        valid(A)[m, n] && return storage(A)[m, n]
    end

    result = integrate(A, o_1, o_2)
    A[o_1, o_2] = result
    return result
end


function Base.setindex!(A::Operator, value, o_1::OnGrid, o_2::OnGrid)
    function add_index(o::OnGrid)
        indices(A)[objectid(o)] = current_index(A)
        current_index!(A, current_index(A) + 1)
    end

    haskey(indices(A), objectid(o_1)) || add_index(o_1)
    haskey(indices(A), objectid(o_2)) || add_index(o_2)

    while size(storage(A), 1) <= current_index(A) 
        expand_storage(A)
    end
    m, n = indices(A)[objectid(o_1)], indices(A)[objectid(o_2)]
    storage(A)[m, n] = value
end

Base.:|(A::Operator, o::OnGrid) = A => o
Base.:|(o::OnGrid, A::Operator) = o => A
Base.:|(o_1::OnGrid, A_and_o_2::Pair{<:Operator, <:OnGrid}) = let (A, o_2) = A_and_o_2
    A[o_1, o_2]
end
Base.:|(o_1_and_A::Pair{<:OnGrid, <:Operator}, o_2::OnGrid) = let (o_1, A) = o_1_and_A
    A[o_1, o_2]
end

# Base.:|(A::Operator, y::AbstractVector) = [A => o_2 for o_2 in y]
Base.:|(x, y::AbstractVector{<:OnGrid}) = [x | o_2 for o_2 in y]
Base.:|(x::AbstractVector{<:OnGrid}, y) = [o_1 | y for o_1 in x]
Base.:|(x::AbstractVector{<:OnGrid}, y::AbstractVector) = [o_1 | o_2 for o_1 in x, o_2 in y]
Base.:|(x::AbstractVector, y::AbstractVector{<:OnGrid}) = [o_1 | o_2 for o_1 in x, o_2 in y]

# Base.:|(A::Operator, x::Vector{OnGrid}, y::Vector{OnGrid}) = [A[o_1, o_2] for o_1 in x, o_2 in y]
# Base.:|(A::Operator, x::Vector{OnGrid}, o_2::OnGrid) = [A[o_1, o_2] for o_1 in x]
# Base.:|(A::Operator, o_1::OnGrid, y::Vector{OnGrid}) = [A[o_1, o_2] for o_2 in y]

function resemble(operator::NumericalOperator{S}, ::Type{T}, new_elements=nothing) where {S <: Grid, T <:Grid}
    g = grid(operator)
    S == dual_grid(T) && (g = transform_grid(g))
    if new_elements === nothing
        new_elements = zeros(eltype(elements(operator)), size(g))
    end
    NumericalOperator(g, new_elements)
end