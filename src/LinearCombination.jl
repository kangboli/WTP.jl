export LinearCombination, basis, basis!, coefficients, coefficients!, add, negate, mul, minus, braket, basis_transform

"""
A linear combination is a set of basis vector and a set of _coefficients with
basic algebra and an inner product (braket) defined.

This abstract type is to provide

1. Operations between linear combinations in different basis set.
2. Simple syntax for these operations.
"""
abstract type LinearCombination end

"""
The basis is always internally stored as kets. 
"""
_basis(linear_combination::LinearCombination) = linear_combination.basis

"""
The basis of a linear comb. Since the basis are stored internally as kets,
they have to be conjugated before returned if the linear comb is a bra.
"""
basis(linear_combination::LinearCombination) = ket(linear_combination) ? _basis(linear_combination) : [dagger(b) for b in _basis(linear_combination)]

_basis!(linear_combination::LinearCombination, new_basis) = linear_combination.basis = new_basis

"""
The _coefficients of the linear combination. 
"""
coefficients(linear_combination::LinearCombination) = linear_combination._coefficients
coefficients!(linear_combination::LinearCombination, new_coefficientsoefficients) = linear_combination.coefficients = new_coefficientsoefficients


"""
Arithmatics. Only minus has a default implementation.
add, mul, negate should be implemented by the concrete types.
"""
minus(l_1, l_2) = add(l_1, negate(l_2))

"""
    braket(bra, ket)

Braket of two linear combinations. The first argument has to be a bra, and the
second must be a ket.  They can be in any basis so long as braket is defined for
all pairs of basis vectors.

⟨ψ|ϕ⟩ = ∑ᵢⱼ aᵢ†bⱼ ⟨αᵢ|βⱼ⟩
"""
function braket(l_1::LinearCombination, l_2::LinearCombination)
    result = 0

    b12 = [braket(b1, b2) for b1 in basis(l_1), b2 in basis(l_2)]
    return transpose(coefficients(l_1)) * b12 * coefficients(l_2)
    # for (b1, c1) in zip(basis(l_1), coefficients(l_1))
    #     for (b2, c2) in zip(basis(l_2), coefficients(l_2))
    #         result +=  c1 * c2 * braket(b1, b2)
    #     end
    # end
    return result
end

Base.:+(l_1::LinearCombination, l_2::LinearCombination) = add(l_1, l_2)
Base.:-(l_1::LinearCombination) = negate(l_1)
Base.:-(l_1::LinearCombination, l_2::LinearCombination) = minus(l_1, l_2)
Base.:*(scalar, l_1::LinearCombination) = mul(scalar, l_1)
Base.:*(vector::AbstractVector, l_1::LinearCombination) = [mul(s, l_1) for s in vector]
Base.:*(l_1::LinearCombination, scalar) = mul(scalar, l_1)
Base.:*(l_1::LinearCombination, vector::AbstractVector) = mul(vector, l_1)
Base.:*(l_1::LinearCombination, l_2::LinearCombination) = braket(l_1, l_2)
Base.adjoint(linear_combination::LinearCombination) = dagger(linear_combination)

"""
Transfrom the linear combination into a different basis.

|ψ⟩ = ∑ᵢ bᵢ|βᵢ⟩ = ∑ aⱼ|αⱼ⟩
a = S⁻¹ R b, Rᵢⱼ = ⟨αᵢ|βⱼ⟩, Sᵢⱼ = ⟨αᵢ|αⱼ⟩
"""
function basis_transform(_coefficients, old_basis, new_basis)
    S = [braket(dagger(ai), aj) for ai in new_basis, aj in new_basis]
    R = [braket(dagger(ai), bj) for ai in new_basis, bj in old_basis]
    return S \ (R * _coefficients)
end


