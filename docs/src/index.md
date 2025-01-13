
```@raw html
<img src="https://github.com/kangboli/WTP.jl/blob/main/wtp_logo.svg?raw=true" width="200"/>
```

## Introduction

`WTP.jl` provides a thin layer of low cost abstractions for processing electronic
orbitals in condense phase.  It is basically a specialized vector bundle library.


## Abstractions

- Periodic and non-orthogonal grids with support for winding. This is a
  unified interface for the Brillouin zone, the reciprocal lattice, the
  homecell, and the crystal lattice.
- Functions defined on these grids (vector bundles) and operation on these
  functions including indexing and FFT. This is an unified interface for 
  the Bloch waves, the Wannier functions, and the $u_{n \mathbf{k}}$ part.


