WTP provides some IO capabilities for interfacing with Quantum Espresso and Wannier90.

This documentation represents my best effort on figuring out the file formats.
It does not serve as a reference, but it's probably as good as any. If there is
a difference from the documentation of Quantum Espresso/Wannier90, you should
rely on their documentation instead (if you can find them :P).


## `wfcX.dat` files

The `wfcX.dat` files are used by Quantum Espresso to store the wave functions in
binary format. It seems that they are moving to HDF5, but I couldn't find any
documentation on that just yet as of 2021.

### File format 

```@docs
WFC
```

```@docs
WFC(::String)
```

### Read a `.save` Directory

```@docs
wave_functions_from_directory
```

```@docs
orbital_set_from_save
```

```@docs
single_orbital_from_wave_functions
```

```@docs
load_evc!
```

### K-point map

```@docs
i_kpoint_map
```

## `UNKXXXXX.1` files

```@docs
UNK
```

```@docs
UNK(::String)
```

```@docs
wannier_from_unk_dir
```

```@docs
single_orbital_from_unk
```

## `.mmn` Files

```@docs
MMN
```

```@docs
MMN(::String)
```

<!-- 
```@docs
NeighborIntegral(::MMN, ::Dict{Int64, KPoint})
``` -->

!!! note 

    For Γ-point calculations, the neighbor integrals are between Γ-point of the first Brillouin zone and that of other (most likely adjacent) Brillouin zones.

## `.amn` Files

```@docs
AMN
```

```@docs
AMN(::String)
```

```@docs
Gauge(::Grid, ::AMN, ::Dict{Int64, KPoint})
```
