A grid vector is a point on our grid. It is primarily used to index functions defined on a grid, but it also behaves like a grid point, which means that you can do arithmatics with them. 

```@setup grid_vector
using WTP
using Gadfly
```

## Creating a Grid Vector

The most typical way of creating a grid vector is by indexing a grid.
Suppose that we have a grid

```@example grid_vector
g₁, g₂, g₃ = CARTESIAN_BASIS
reciprocal_basis = (g₁, 2g₂, 3g₃)
sizes = (4, 4, 4)
lattice = ReciprocalLattice3D(reciprocal_basis, size_to_domain(sizes))
```

We can obtain a grid vector by indexing into it

```@example grid_vector
g = lattice[0, 1, 2]
```

This is equivalent to an explicit construction

```@example grid_vector
g = GridVector{ReciprocalLattice}(lattice, [0, 1, 2], true)
```

## Indexing a Grid Vector

The primary usage of the grid vector is to index functions defined on our `Grid`.
Consider some planewave defined on the reciprocal lattice (a $\delta$ function)

```@example grid_vector
ψ = map(g->(g==lattice[0,0,0]) |> Float64, lattice)
ψ[lattice[1, 0, 0]], ψ[lattice[0, 1, 0]]
```

!!! note 

    If your function is on a different grid, you would need a different type of grid vector. This restriction is intentional to prevent abuses of indices.


This type of indexing make code far more readable than it's Fortran
counterparts.

Indexing with ranges is also supported
```@example grid_vector
spy(ψ[lattice[:, :, 0]])
ans |> SVG("planewave.svg", 4inch, 4inch); nothing # hide
```

![](planewave.svg)

It is often time desirable to treat a function as a raw 1D array so that matrix
methods apply. This is not hard. `elements(ψ)` can be treated as a 1D array.
```@example grid_vector
elements(ψ)[1:3]
```

The harder part is to find out which grid vectors these indices correspond to.
This conversion is made simple by indexing the grid with a different syntax.
```@example grid_vector
ψ[lattice(1:3)]
```

Converting from a grid vector to an 1D index is discouraged since it is a bad idea most of the time, but it can be done if you need it.
```@example grid_vector
i = lattice[1, 0, 0] |> 
(g -> miller_to_standard(g, [0, 0, 0])) |>
(v -> three_to_one(v..., size(lattice)))
```
The vector `[0, 0, 0]` is for a feature that is not yet implemented.


## Arithmatics 

Arithmatics work mostly as one expects. This makes it very easy to 
find neighbors and traverse a grid.

### Add (and Substract)

```@example grid_vector
lattice[1, 0, 0] + lattice[0, 0, 1]
```

### Scalar Multiply

```@example grid_vector
2lattice[1, 0, 0] 
```

### Inner Product

```@example grid_vector
homecell = transform_grid(lattice)
homecell[1, 0, 0]' * lattice[1, 0, 0]
```

### Cartesian Coordinates

```@example grid_vector
cartesian(homecell[1, 0, 1])
```

### Norm

```@example grid_vector
using LinearAlgebra
norm(homecell[1, 0, 1])
```

## Overflow 

Grid vectors can reach out of the domain of the underlying grid, and this is allowed.

```@example grid_vector
r = homecell[3, 7, -9]
```

You can check if it is out of the grid 

```@example grid_vector
has_overflow(r)
```

and by how much (in units of entire grids)

```@example grid_vector
overflow(r)
```

To find its periodic image within the grid

```@example grid_vector
reset_overflow(r)
```