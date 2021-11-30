WTP provides some IO capabilities for interfacing with Quantum Espresso and 
Wannier90

```@setup iopw
using WTP
```

## `wfcX.dat` files

The `wfcX.dat` files are used by Quantum Espresso to store the wave functions in
binary format. It seems that they are moving to HDF5, but I couldn't find any
documentation on that just yet as of 2021.

### File format 

```@docs
WFC
```

### Example 

For an Si example, suppose that you are running Julia in the project root of `WTP`, some data is available for testing in the following dir.

```julia
path_to_si = "../test/test_data/test_1"
```

To read a wfc file, just pass the path to the file to `WFC`.

```@example iopw
path_to_si = "../../test/test_data/test_1" # hide
wave_function = WFC(joinpath(path_to_si, "si.save/wfc1.dat"))
wave_function.n_band, wave_function.gamma_only, wave_function.n_planewaves
```

### Read a `.save` Directory

The `wfcX.dat` files only gives you the content of a single wave function, but what you really want is probably reading the entire `.save` directory and convert them to orbitals on grids.
To do that,

```@example iopw
# Read the entire directory.
wave_functions_list = wave_functions_from_directory(joinpath(path_to_si, "si.save"))
# Convert them to orbitals.
u = wannier_from_save(wave_functions_list)
typeof(u)
```

The result is a function on a `BrillouinZone3D`, where the value at each grid point
is a list of orbitals.

```@example iopw
Γ = grid(u)[0, 0, 0]
typeof(u[Γ])
```

## `UNKXXXXX.1` files

### File format 

The `UNKXXXXX.1` files are produced by `pw2wannier90.x` to show orbitals in real space.

```@docs
UNK
```


### Example

To read an `UNKXXXXX.1` file, just pass in the file path to `UNK`.

```@example iopw
path_to_si_save = "../../test/test_data/test_1/si.save" # hide
wave_function = UNK(joinpath(path_to_si, "unk/UNK00001.1"))
wave_function.n_band, wave_function.nx
```

### Read an `unk` Directory

Similar to the `wfcX.dat` files, you probably want to load an entire directory.
Since the `UNKXXXXX.1` does not carry much information, we need the `wave_functions_list` from 
the `.save` directory as well.

```@example iopw
# u_real = wannier_from_unk_dir(joinpath(path_to_si, "unk"), wave_functions_list)
# typeof(u_real)
```

The result is identical from reading from the `wfcX.dat` files, and we will not show it here.

## `.mmn` Files

The `.mmn` files are used by Wannier90 for storing the pairwise integrals of the
orbitals correponding to neighboring k-points. To read a `.mmn` file, pass in
the file path to `MMN`.

```@example iopw
mmn = MMN(joinpath(path_to_si, "output/pw2wan/si.mmn"))
mmn.n_kpoint, mmn.n_neighbor
```

We like to convert this raw data to `NeighborIntegral` before looking up
integrals from it. One problem that we encounter is that the k-points are stored
as a single integer in the `.mmn` files. To find out which k-points these
integers correspond to, we have to construct a mapping using the `wfcX.dat`
files.

```@example iopw
k_map, brillouin_zone = i_kpoint_map(wave_functions_list)
k_map[1], k_map[2], k_map[3]
```

This is enough information for us to build a `NeighborIntegral`. Indexing 
it with a pair of k-points gives a pairwise integral matrix between the two k-points.

```@example iopw
neighbor_integral = NeighborIntegral(mmn, k_map)
neighbor_integral[Γ, brillouin_zone[1, 0, 0]][1:3, 1:3]
```

!!! note 

    For Γ-point calculations, the neighbor integrals are between Γ-point of the first Brillouin zone and that of other (most likely adjacent) Brillouin zones.



## `.amn` Files

The `.amn` files are used by Wannier90 for storing the gauge transform.
Reading it is similar to reading a `.mmn` file.

```@example iopw
amn = AMN(joinpath(path_to_si, "output/pw2wan/si.amn"))
amn.n_kpoint
```

An `AMN` object also just stores the raw data, and we need to convert this to a
`Gauge` to use it. Like `.mmn` files, we need a mapping of k-points as well as the Brillouin zone to construct the gauge, which can then be indexed with k-points.

```@example iopw
U = Gauge(brillouin_zone, amn, k_map)
U[Γ][1:3, 1:3]
```

