# Orbital set

```@docs
OrbitalSet
```


## Create an orbital set

```@docs
init_orbital_set(::BrillouinZone, ::Type{T}) where T<:OnGrid
```

## Basic Usage
### Gauge

```@docs
Gauge
```

```@docs
Gauge(::T, ::Integer) where T <: BrillouinZone
```

```@docs
gauge
```

```@docs
set_gauge
```

```@docs
commit_gauge
```

### Dimensions

```@docs
orbital_grid
```

```@docs
superlattice
```

```@docs
supercell
```

```@docs
n_band(::OrbitalSet)
```

### FFT

```@docs
fft(::OrbitalSet{UnkBasisOrbital{T}}) where {T<:HomeCell}
```

```@docs
ifft(::OrbitalSet{UnkBasisOrbital{T}}) where {T<:ReciprocalLattice}
```

In place fft is also available through `@ifft` and `@fft`.
The argument to these macros will be set to nothing.

## Indexing

### Single indexing

```@docs
getindex(::OrbitalSet, ::KPoint)
```

### Array Indexing

Not yet implemented.

### Double indexing

```@docs
getindex(::OrbitalSet, ::Integer, ::KPoint)
```

### Callable indexing

```@docs
wannier_orbitals
```

```@docs
reciprocal_densities
```
