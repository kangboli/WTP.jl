The orbitals are defined somewhat abstractly as functions on grids. The grid
points are often times either real space points or wave numbers, but more generally they are sets of parameters of the basis functions.

```@setup orbital
using WTP
```

## Creating an Orbital

There are two ways of creating an orbital. The primary way is to read from the output of Quantum Espresso, but that goes into the IO part of the package. 
The other way is to map a function onto a grid.

Suppose that we have a reciprocal grid
```@example orbital
g₁, g₂, g₃ = CARTESIAN_BASIS
reciprocal_basis = (g₁, 2g₂, 3g₃)
sizes = (4, 4, 4)
lattice = ReciprocalLattice3D(reciprocal_basis, size_to_domain(sizes))
```

We can create a single plane wave as
```@example orbital
gₓ = lattice[1,0,0]
ψ = map(g->(g==gₓ) |> Float64, lattice)
ψ[gₓ]
```

We could also do this in the home cell.

```@example orbital
homecell = transform_grid(lattice)
ϕ = map(r->exp(1im * gₓ' * r), homecell) |> wtp_normalize!
(x->round(x, digits=3)).(ϕ[homecell[-2:1,0,0]])
```

## Fast Fourier Transform

The package is designed to make this as easy as it can get. To perform a forward FFT, 
```@example orbital
ϕ̂ = fft(ϕ)
ϕ̂[gₓ]
```
The type of the orbital and the underlying grid are also transformed accordingly.
```@example orbital
typeof(ϕ), typeof(ϕ̂)
```

The inverse FFT is,
```@example orbital
ψ̂ = ifft(ψ)
ψ̂[homecell[-2:1, 0, 0]]
```

## Indexing an Orbital

See the section on grid vectors.

## Arithmatics

### Addition and Subtraction

```@example orbital
using LinearAlgebra
norm(elements(ψ + ϕ̂)), norm(elements(ψ - ϕ̂))
```

### Inner Product

This is a very common operation, so the syntax must be simple.

```@example orbital
ψ' * ϕ̂
```

If you like to custom the inner product, the method to overload is
`braket`.

```@example orbital
braket(ψ', ϕ̂)
```

### Translation

A translation is useful for efficient phase shifting.
To translate an orbital by $[-1, 0, 0]$ (move the orbital to the left by one grid point), we have

```@example orbital
(ψ >> [-1, 0, 0])[lattice[0, 0, 0]]
```

## Linear Combination of Orbitals (Symbolic)

Given a few basis orbitals, we can construct a linear combination of them for a basic set of symbolic manipulation.  A linear combination 
is named `UnkOrbital` due to my stupidity.

```@example orbital
planewave(wave_numbers...) = map(g->(g==lattice[wave_numbers...]) |> 
Float64, lattice) |> UnkOrbital

ϕ₁ = planewave(1, 0, 0)
ϕ₂ = planewave(0, 1, 0)
ϕ₃ = planewave(0, 0, 1)
ψ = 1ϕ₁ + 2ϕ₂ + 3ϕ₃
```

Arithmatics should work as one expects

### Addition and Subtraction

```@example orbital
ψ - 2ϕ₁ - ϕ₃
```

### Inner Product

```@example orbital
(ψ - 2ϕ₁ - ϕ₃)' * ψ 
```

### Matrix Multiplication

```@example orbital
[ϕ₁, ϕ₂, ϕ₃]' * [0 0 1;
                 0 1 0;
                 1 0 0] * [ϕ₁, ϕ₂, ϕ₃]

```
