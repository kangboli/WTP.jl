export ket, ket!, dagger, braket, Vector3, CARTESIAN_BASIS,
vector3_to_matrix, matrix_to_vector3

"""
A Basis vector is an object that supports inner product with other basis
vectors.  Basis is not defined as a type because of the lack of multiple
subtyping.  One just has to guess if something is a basis.
"""

"""
    ket(b)

Returns whether the basis vector is a ket.
A basis can be either a bra or a ket.
"""
ket(b) = b.ket
ket!(b, is_ket::Bool) = b.ket = is_ket

"""
The conjugate of a basis vector.
"""
dagger(b) = @set b.ket = !ket(b)    
function dagger!(b) 
    b = ket!(b, !ket(b))
end

"""
The braket of a bra and a ket should be a number.
"""
function braket(b_1, b_2)::Number
    !ket(b_1) && ket(b_2) || error("braket should have a bra and a ket.")
    return b_1.vector' * b_2.vector
end


"""
A wrappper over 3D static vectors as basis vectors.
"""
struct Vector3{T <: Number}
    vector::SVector{3, T}
    ket::Bool
end

Vector3(x, y, z, ket=true) = Vector3(SVector(x, y, z), ket)

"""
Convert a set of Vector3 to a matrix, whose columns are the vectors.
"""
vector3_to_matrix(v3s) = hcat((v3 -> v3.vector).(v3s)...)

"""
Convert a matrix to a set of Vector3, which are the columns of the matrix.
"""
matrix_to_vector3(mat, ket=true) = Tuple(Vector3(col..., ket) for col in eachcol(mat))

function Base.show(io::IO, v::Vector3) 
    ket(v) ? @printf(io, "ket: ") : @printf(io, "bra: ")
    @printf(io, "%.3f, %.3f, %.3f", v.vector...)
end

function Base.:(==)(v_1::Vector3, v_2::Vector3)
    return isapprox(v_1.vector, v_2.vector, atol=1e-7) && ket(v_1) == ket(v_2)
end

const CARTESIAN_BASIS = (Vector3(1., 0., 0.), Vector3(0., 1., 0.), Vector3(0., 0., 1.))

Base.:*(s::Number, v::Vector3) = Vector3(s*v.vector, ket(v))
Base.:*(v::Vector3, s::Number) = Vector3(s*v.vector, ket(v))
Base.:/(v::Vector3, s::Number) = Vector3(v.vector/s, ket(v))
