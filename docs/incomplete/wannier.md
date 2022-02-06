The wannier functions are somewhat complicated to represent, since it does not
fit into any of our grid, but it can be represented as a linear combination of
ones that do fit, which are the $u_{nk}(r)$ orbitals. We will provide a
structure `Wannier` for organizing and manipulating these orbitals.

## Setup 

We will read in some test data as we did in one of the sections of File IO.

```julia
path_to_si = "../test/test_data/test_1"

using WTP
wave_functions_list = wave_functions_from_directory(joinpath(path_to_si, "si.save"))
u = wannier_from_save(wave_functions_list)
k_map, brillouin_zone = i_kpoint_map(wave_functions_list)
```

```@setup wannier
using WTP

path_to_si = "../../test/test_data/test_1" # hide
wave_functions_list = wave_functions_from_directory(joinpath(path_to_si, "si.save"))
u = wannier_from_save(wave_functions_list)
k_map, brillouin_zone = i_kpoint_map(wave_functions_list)
Γ = brillouin_zone[0, 0, 0]
```

## Gauge

The structure `Wannier` contains a `Gauge`, which we will discuss first.

```@docs
Gauge
```

There are two ways of constructing a gauge.  The first way is to construct an
empty gauge from a Brillouin zone, and set the matrices at each k-point
manually.  Here, empty means that the matrices are `undef`, not identity
matrices.

```@example wannier
U = Gauge(brillouin_zone)
isdefined(elements(U), 1)
```

Iterating over all the k-points is straight-forward.

```@example wannier
using LinearAlgebra

for k in brillouin_zone
    U[k] = Diagonal(ones(ComplexF64, 3))
end

U[Γ]
```

The other way of constructing a gauge is from an `AMN` object.

```@example wannier
amn = AMN(joinpath(path_to_si, "output/pw2wan/si.amn"))
U = Gauge(brillouin_zone, amn, k_map)
size(U[Γ])
```

## Wannier

```@docs
Wannier
```

It is not straightforward to construct a `Wannier` from scrach since we would have to make up all of the $u_{nk}$ orbitals, but constructing it from a `.save` directory is easy, and we already did it.

```julia
wave_functions_list = wave_functions_from_directory(joinpath(path_to_si, "si.save"))
u = wannier_from_save(wave_functions_list)
```

Since `Wannier` is a function on a grid (the Brillouin zone), it can be indexed with a `GridVector` $k$.
The result is the set of $u_{nk}$ orbitals at $k$.

```@example wannier
typeof(u[Γ])
```

However, to get the $u_{nk}$ orbital together with the gauge transform, 
the syntax is `u[n, k]`.

```@example wannier
typeof(u[1, Γ])
```

It is tempting to construct the Wannier functions explicitly, but it turns out
that doing so is mostly useless, although possible.

```@example wannier
sum(brillouin_zone) do k
    (u[1, k] |> UnkBasisOrbital) >> coefficients(k) |> UnkOrbital
end
```

This code is a bit complicated. Since `u[1, k]` is a `UnkOrbital` (linear
combination), we need to consolidate it to a `UnkBasisOrbital` before we can
apply as phase shift.

## Neighbor Integrals

The neighbor integrals roughly corresponds to the $M_{mn}^{k, b}$ matrices
in Wannier90, but is implemented with a dictionary hashing on pairs of 
k-points (`GridVector{<:BrillouinZone}`) for more flexibility and readibility.

There are generally two ways to construct a neighbor integral. The primary way
is to construct a finite difference scheme, which we will discuss later in
detail. The other way is to read from a `.mmn` file. 

```@example wannier
mmn = MMN(joinpath(path_to_si, "output/pw2wan/si.mmn"))
M = NeighborIntegral(mmn, k_map)
M[Γ, brillouin_zone[1, 0, 0]][1, 1]
```

We can apply a gauge transform to a neighbor integral directly, so that 
we don't have to recompute the integrals every time we transform the orbitals.

```@example wannier
M̃ = gauge_transform(M, U)
M̃[Γ, brillouin_zone[1, 0, 0]][1, 1]
```

## Finite Difference

WTP currently implements the exact same finite difference scheme as Wannier90 to
reproduce the same center and spread. To create a scheme, we only need the `Wannier`
object. 

```@example wannier
scheme = W90FiniteDifference3D(u)
shells(scheme), weights(scheme)
```

This will find shells of nearest neighbors, compute their weights, and construct
a `NeighborIntegral` (takes some time). You can also find the shells and the weights separtely if you want them without having to do the integrals.

### Find the Shells

Finite difference requires finding shells of nearest neighbors of a k-point $k$. 
The neighbors in each shell has the same distance to $k$. We provide a function for finding shells of nearest neighbors for the Γ-point.

```julia
shells = find_shells(brillouin_zone, 2)
```

### Find the Weights 

Each shell in the Finite difference may have a different weight.
Finding the weights requires solving a linear equation. The number of shells 
required varies by the problems. In our case, we need only one.

```julia
weights = compute_weights(shells[1:1])
```

## Center and Spread

To compute the center and the spread of a Wannier function, we need 

1. The neighbor integral transformed by the gauge.
2. The finite difference scheme.
3. The band index of the Wannier function.

The center (bohr) 

```@example wannier
c = center(M, scheme, 1)
```

The spread (bohr²)

```@example wannier
σ² = second_moment(M, scheme, 1) - norm(c)^2
```
