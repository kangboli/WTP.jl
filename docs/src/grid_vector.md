
```docs
AbstractGridVector
```

## Creating a grid vector

```@docs
make_grid_vector
```

## Basic usage

### The underlying grid.

```@docs
grid(::AbstractGridVector)
```

```@docs
set_grid
```

### Basis and coefficients

```@docs
basis(::AbstractGridVector)
```

```@docs
coefficients(::AbstractGridVector)
```

## Overflow

```@docs
overflow
```

```@docs
has_overflow
```

```@docs
wrapped
```

```@docs
reset_overflow
```

```@docs
add_overflow
```

## Arithmatic

```@docs
add(::T, ::T) where T <: AbstractGridVector
```

```@docs
negate(::AbstractGridVector)
```

```@docs
mul(::Int, ::AbstractGridVector)
```

```@docs
coordinates(::AbstractGridVector)
```

```@docs
LinearAlgebra.norm(::AbstractGridVector)
```

```@docs
braket_
```
