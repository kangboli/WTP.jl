Γ-point calculations had their own branch of `if` statements in Wannier90,
but it turns out that this wasn't necessary, not even for performance reasons. We will use Benzene as an example

## Reading from `.save`

```@setup gamma
using WTP
const path_to_benzene = "../../test/test_data/test_3"
```

```julia
const path_to_benzene = "../test/test_data/test_3"
```

```@example gamma
wave_functions_list = wave_functions_from_directory(joinpath(path_to_benzene, "benzene.save"))
u = orbital_set_from_save(wave_functions_list);
k_map, _ = i_kpoint_map(wave_functions_list)
grid(u)
```

## Finite Difference

For Γ-point calculations, the center and the spread are evaluated with a
separate procedure that is not evidently finite difference. However,
the dedicated procedure is exactly equivalent to a finite difference
scheme if we consider the neighbors to be the Γ-point of adjacent Brillouin zones.

```@example gamma
amn = AMN(joinpath(path_to_benzene, "output/pw2wan/benzene.amn"))
U = Gauge(grid(u), amn, k_map)
scheme = CosScheme3D(u, 1)
M = gauge_transform(neighbor_basis_integral(scheme), U)

center(M, scheme,  1)
```