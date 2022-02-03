WTP provides *non-orthogonal 3D periodic* grids that support

1. Three Indexing schemes and iteration.
2. Fourier transform.

The goals is to provide a uniform interface for dealing with

1. The reciprocal lattice.
2. The Brillouin zone.
3. The crystal lattice.
4. The homecell.

```@setup grid
using WTP
using Gadfly
```
## Centering Convention

When dealing with multiple periodic grids, two questions on the matter of
convention give people a lot of headache. Consider the case of 1D,
the questions are

1. What are the left-most and the right-most grid points?
2. Which grid point is considered the center?

We will deal with these problems by imposing a simple convention for all grids in our package, and we refer to it as the centering convention. Again, in the case of 1D, suppose that the step size of the grid is 
$a$.

- The grid is even.
- For a grid with $2 N$ grid points, the leftmost grid point is at $-N a$, and the rightmost grid point is at $(N-1) a$.
- The center is at the origin.

## Creating a Grid

To make a home cell with basis vectors

$$a_1 = \begin{pmatrix}
0\\ \sqrt{3}/2 \\ 0
\end{pmatrix} \quad
a_2 = \begin{pmatrix}
1\\ 1/2 \\ \sqrt{3}/2 
\end{pmatrix} \quad
a_3 = \begin{pmatrix}
0\\ 0 \\ 1/2
\end{pmatrix}$$

and domain $(-10, 9)$ in each direction in terms of the basis vectors.

```@example grid
homecell_basis = (Vector3(0, sqrt(3)/2, 0), Vector3(1, 1/2, sqrt(3)/2), Vector3(0, 0, 1/2))
homecell = make_grid(HomeCell3D, homecell_basis, ((-10, 9), (-10, 9), (-10, 9)))
```

We provide four different grids by default

```@docs
HomeCell3D
ReciprocalLattice3D
BrillouinZone3D
RealLattice3D
```

## Indexing a Grid

A grid is primarily indexed by its Miller indices in a general sense. The result of the 
indexing is a `GridVector` parametrized with the grid itself.

```@example grid
origin = homecell[0, 0, 0]
```

Ranges are supported, and the result will be a multi-dimensional array.
```@example grid
homecell[-1:0, 0, -1:0]
```

The grid can also be indexed with a single number using parenthesis.  This more
or less vectorizes the grid for matrix algorithms. There are routines for
converting this number to its equivalent 3D index and vice versa.

```@example grid
homecell(4)
```

```@example grid
indices_3d = standard_to_miller(size(homecell), one_to_three(4, size(homecell)), [0, 0, 0])
homecell[indices_3d...]
```

!!! warning "Do not pass around an integer as 3D indices"
    
    One should not pass an integer as a 3D index, amongst many other things
    people abuse integers for.  An integer does not have the grid as its context, so
    one should pass a `AbstractGridVector` instead most of the time.

## Iterations and Maps over Grids

All grids are iterable. For example, if we like to sum over the norm of every
grid point in the homecell,

```@example grid
sum(homecell) do r
    r' * r
end
```

We can also map over a grid, the result will be a `SimpleFunctionOnGrid`.
For example, we can generate a Gaussian on the home cell, and visualize a slice of it.

```@example grid
gaussian = map(r->exp(-r'*r), homecell)
spy(gaussian[homecell[0, :, :]])
ans |> SVG("gaussian.svg", 4inch, 4inch); nothing # hide
```

![](gaussian.svg)

The result does not look isotropic even though our Gaussian is. This is because the grid is not orthogonal. Keep this in mind for all the pictures.

## Transform a Grid

If we apply a Fourier transform on a function (see the section on orbitals)
defined on a grid, the result is a function defined on a different grid, which
we refer to as the dual grid.

```@example grid
guassian_reciprocal = fft(gaussian)
reciprocal_lattice = grid(guassian_reciprocal)
spy(abs.(guassian_reciprocal[reciprocal_lattice[0, :, :]]))
ans |> SVG("gaussian_reciprocal.svg", 4inch, 4inch); nothing # hide
```
![](gaussian_reciprocal.svg)

```@example grid
grid(guassian_reciprocal)
```

One can also explicitly transform only the grid.
```julia
reciprocal_lattice = transform_grid(homecell)
```

`HomeCell3D` and `ReciprocalLattice3D` are dual grids and
`BrillouinZone3D` and `RealLattice3D` are dual grids.
## Some other Methods

```@docs
domain(::Grid)
array(::Grid)
n_dims(::Type{<:Grid})
```