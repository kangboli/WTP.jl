using Coverage, Dates

export indent, size_to_domain, miller_to_standard, standard_to_miller, three_to_one, one_to_three, @set, process_coverage

"""
    size_to_domain(sizes)

Convert a set of three sizes into a centered domain following the centering
convention.

Example: 

```@jldoctest miller
julia> size_to_domain((3, 4, 5))
((-1, 1), (-2, 1), (-2, 2))
```
"""
function size_to_domain(sizes)
    domain(s::Integer) = iseven(s) ? (-s ÷ 2, s ÷ 2 - 1) : (-s ÷ 2, s ÷ 2)

    return Tuple(domain(n) for n in Int.(sizes))
end

function indent(str::String)
    lines = split(str, "\n")
    lines = [join(fill(" ", 4)) * line for line in lines]
    return join(lines, "\n")
end


# """
#     miller_to_standard(sizes, indices)

# Convert a set of miller indices to standard indices without offset. 
# """
# miller_to_standard(
#     sizes::Tuple,
#     indices::Tuple,
# ) = collect(
#     Iterators.map(
#         (m, s) -> m >= 0 ? m + 1 : m + s + 1,
#         indices,
#         sizes,
#     ),
# )

"""
    miller_to_standard(sizes, indices, offsets)

Convert a set of miller indices to standard indices. 

Example:

```@jldoctest miller
julia> miller_to_standard((4, 4, 4), (0, 0, 0), (0, 0, 0))
(1, 1, 1)
julia> miller_to_standard(sizes, (-3, -1, 1), (1, 0, 0))
(3, 4, 2)
```
"""
miller_to_standard(sizes::Tuple, indices::Tuple, offsets::Tuple) = 
    map(miller_to_standard_1d, indices, offsets, sizes)

miller_to_standard_1d(m::Int, o::Int, s::Int) = m + o >= 0 ? m + o + 1 : m + o + s + 1


# """
#     standard_to_miller(sizes, indices)

# Convert a set of standard indices to miller indices without offsets.
# """
# standard_to_miller(
#     sizes::Tuple,
#     indices::Tuple,
# ) = collect(
#     Iterators.map(
#         (m, s) -> let z = iseven(s) ? s ÷ 2 : s ÷ 2 + 1
#             m <= z ? m - 1 : m - 1 - s
#         end,
#         indices,
#         offsets,
#         sizes,
#     ),
# )

"""
    standard_to_miller(sizes, indices, offsets)

Convert a set of standard indices to miller indices.

Example:

```@jldoctest miller
julia> standard_to_miller(sizes, (1, 1, 1), (0, 0, 0))
(0, 0, 0)
julia> standard_to_miller(sizes, (3, 4, 2), (1, 0, 0))
(-3, -1, 1)
```
"""
standard_to_miller(sizes::Tuple, indices::Tuple, offsets::Tuple) = 
    map(standard_to_miller_1d, indices, offsets, sizes)

standard_to_miller_1d(m::Int, o::Int, s::Int) = let z = iseven(s) ? s ÷ 2 : s ÷ 2 + 1
    return m <= z ? m - 1 - o : m - 1 - s - o
end


"""
    one_to_three(i, sizes)

Convert 1D indices to 3D ones.

The conversion is best characterized as the following relation.

```julia
reshape(A, prod(sizes))[i] == A[one_to_three(i, sizes)...]
```

That is, indexing with a 3D index and a 1D index should be identical 
if the two indices are mapped to each other with `one_to_three`.

Example:

```@jldoctest miller
julia> one_to_three(5, (4, 4, 4))
(1, 2, 1)
```
"""
one_to_three(i::Integer, sizes::NTuple{3, Int}) =
    let (nx, ny, _) = sizes
        i = i - 1
        iz, yx = i ÷ (ny * nx) + 1, rem(i, ny * nx)
        iy, ix = yx ÷ nx + 1, rem(yx, nx) + 1
        (ix, iy, iz)
    end

"""
    three_to_one(x, y, z, sizes)

Convert 3D indices to 1D ones.

The conversion is best characterized as the following relation.

```julia
reshape(A, prod(sizes))[three_to_one(i, j, k, sizes)] = A[i, j, k]
```

Example: 

```@jldoctest miller
julia> three_to_one(1, 2, 1, (4, 4, 4))
5
```
"""
three_to_one(x::Integer, y::Integer, z::Integer, sizes::NTuple{3, Int}) =
    let (nx, ny, _) = sizes
        (z - 1) * ny * nx + (y - 1) * nx + x
    end


"""
    @set a.b = c

Create a copy of `a` with its field `b` set to `c`. 

This is a castrated version of `@set` from `Setfield`.
"""
macro set(assignement)
    target = assignement.args[1]
    object = target.args[1]
    modify_field = target.args[2]
    value = assignement.args[2]

    esc(:(set($object, $modify_field, $value)))
end

function set(object, modify_field::Symbol, value)
    values = [f == modify_field ? value : getfield(object, f) for f in fieldnames(typeof(object))]
    return typeof(object)(values...)
end

function process_coverage()
    coverage = process_folder("src")
    covered, total = get_summary(coverage)
    @printf("Coverage is at %.2f percent\n", 100*covered/total)
    LCOV.writefile("coverage/lcov_$(now()).info", coverage)
    clean_folder("src")
    clean_folder("test")
end