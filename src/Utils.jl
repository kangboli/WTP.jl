using Coverage, Dates

export indent, size_to_domain, miller_to_standard, standard_to_miller, three_to_one, one_to_three, @set, process_coverage

"""
    size_to_domain(sizes)

Convert a set of three sizes into a domain.

sizes must be an iterable of numbers.
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

"""
## The Miller Indices.

TODO: Offsets have not been tested. 

Convert indices from miller indices to storage indices.  for example, for an
even grid:

miller:   -3 -2 -1  0  1  2  
standard:  4  5  6  1  2  3  

With an offset of 1, we effectively translate the orbital to the right by one.

miller:    -4 -3 -2 -1  0  1 
standard:   5  6  1  2  3  4

For an odd grid, without an offset.

miller:   -3 -2 -1  0  1  2  3
standard:  5  6  7  1  2  3  4

With an offset, we have

miller:   -4 -3 -2 -1  0  1  2
standard:  5  6  7  1  2  3  4

Convention:

Without any offset, for an even grid, the miller indices go from -n+1 to n.
For an odd grid, the miller indices go from -n to n.

To have a different range of miller indices. add offset o. For an even grid,
the indices with offsets go from -n+1-o to n-o. For an odd grid, the indices
with offsets go from -n-o to n-o.
"""


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

>>> reshape(A, prod(sizes))[i] = A[one_to_three(i, sizes)...]

That is indexing with a 3D index and a 1D index should be identical 
if the two indices are mapped to each other with the following functions.
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

>>> reshape(A, prod(sizes))[three_to_one(i, j, k, sizes)] = A[i, j, k]
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