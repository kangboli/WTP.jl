# Grid

```@docs
Grid
```

## Centering convention

When dealing with multiple periodic grids, two questions on the matter of
convention give people a lot of headache. Consider the case of 1D,
the questions are

1. What are the left-most and the right-most grid points?
2. Which grid point is considered the center?

We will deal with these problems by imposing a simple convention for all grids in our package, and we refer to it as the centering convention. Again, in the case of 1D, suppose that the step size of the grid is 
$a$.

- For a grid with $2 N$ grid points, the leftmost grid point is at $-N a$, and the rightmost grid point is at $(N-1) a$.
- For a grid with $2 N + 1$ grid points, the leftmost grid point is at $-N a$, and the rightmost grid point is at $N a$.
- The center is at the origin.

A grid that does not comply to the centering convention is said to have an
offset, which is the same as the center of the grid (offset from the origin).

## Creating a grid

```@docs
make_grid
```

We provide four different grids by default

```@docs
HomeCell3D
ReciprocalLattice3D
BrillouinZone3D
RealLattice3D
```

## Basic usage

### Basis

```@docs
basis(::Grid)
```

```@docs
basis_matrix
```

### Domain

```@docs 
domain
```

```@docs
domain_matrix
```

```@docs
set_domain
```

```@docs
expand(::Grid)
```

### Boundaries

```@docs
mins
```

```@docs
maxes
```

### Center 

```@docs
center(::Grid)
```

### Dimension

```@docs
n_dims
```

```@docs
size(::Grid)
```

```@docs
length(::Grid)
```


### Snap

```@docs 
snap
```

## Indexing a grid

### Single index

```@docs 
getindex(::Grid, ::Integer, ::Integer, ::Integer)
```

### Range index

```@docs
getindex(::Grid, ::AbstractVector{<:Integer}, ::Integer, ::AbstractVector{<:Integer})
```

### Linear index

```@docs
linear_index(::Grid, ::Integer)
```

!!! warning "Do not pass around the 1D indices"
    
    The correspondance between the one dimensional index and the grid point is
    an implementation detail, and one must not rely on it. In particular, one
    should not use it as a function argument.  The two way mapping will be lost
    since indices by themselves have no context (the grid). one should pass a
    `AbstractGridVector` instead most of the time. 

    
## Iterations and maps over grids

```@docs 
iterate(::Grid)
```

## Transform a grid

```@docs
transform_grid
```


