# $U_{n, \mathbf{k}}$ Orbitals

```@docs
AbstractUnkOrbital
```

## Basis orbitals

A `UnkBasisOrbital` inherits `OnGrid` with a few additional properties.

```@docs
UnkBasisOrbital(::T, ::AbstractArray, ::AbstractGridVector{<:BrillouinZone}, ::Integer) where T 
```

```@docs
kpoint(::UnkBasisOrbital)
```

```@docs
kpoint!(::UnkBasisOrbital::UnkBasisOrbital, ::KPoint)
```

```@docs
index_band(::UnkBasisOrbital)
```

```@docs
dagger(::UnkBasisOrbital)
```

```@docs
dagger!(::UnkBasisOrbital)
```

## Linear combination of orbitals.

```@docs
UnkOrbital
```

```@docs
UnkOrbital(::Any)
```

```@docs
dagger(::UnkOrbital)
```

```@docs
dagger!(::UnkOrbital)
```

```@docs
index_band(::UnkOrbital)
```

```@docs
orthonormal(::UnkOrbital)
```

## Indexing of linear combinations

Not yet implemented.

## Semi-symbolic arithmatics (experimental)

```@docs
add(::UnkOrbital, ::UnkOrbital)
```

```@docs
negate(::UnkOrbital)
```

```@docs
mul(::Number, ::UnkOrbital)
```

```@docs
braket(::UnkOrbital, ::UnkOrbital)
```

```@docs
zeros(::UnkOrbital, Vararg)
```