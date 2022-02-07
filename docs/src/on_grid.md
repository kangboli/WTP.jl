```@docs
OnGrid
```

## Creating an Orbital

```@docs
map(::Function, ::T) where T <: Grid
```

## Indexing

### Indexing with a single grid vector

```@docs
getindex(::OnGrid, ::AbstractGridVector)
```

```@docs
setindex!(::OnGrid, ::Any, ::AbstractGridVector)
```

### Indexing with an array of grid vectors

```@docs
getindex(::OnGrid, ::AbstractArray{<:AbstractGridVector})
```

```@docs
setindex!(::OnGrid, ::AbstractArray, ::AbstractArray{<:AbstractGridVector})
```

## Normalization

```@docs
norm(::OnGrid)
```

```@docs
square_normalize
```

```@docs
square_normalize!
```

## FFT

```@docs
fft(::OnGrid{T}) where T <: HomeCell
```

```@docs
ifft(::OnGrid{T}) where T <: ReciprocalLattice
```

```@docs
@fft!
```

```@docs
@ifft!
```


## Arithmatics

```@docs
add(::OnGrid{T}, ::OnGrid{T}) where {T<:Grid}
```

```@docs
negate(::OnGrid{T}) where {T<:Grid}
```

```@docs
minus(::OnGrid{T}, ::OnGrid{T}) where {T<:Grid}
```

```@docs
mul(::OnGrid{T}, ::OnGrid{T}) where {T<:Grid}
```

```@docs
abs2(::OnGrid{T}) where {T<:Grid}
```

```@docs
braket(::OnGrid, ::OnGrid)
```


## Density center

```@docs
compute_r2
```

```@docs
center_spread
```


## Misc

```@docs
expand(::OnGrid)
```

```@docs
vectorize(::OnGrid)
```

```@docs
sparsify
```

```@docs
Base.:(>>)(::OnGrid{S}, ::AbstractGridVector{S}) where {S}
```
