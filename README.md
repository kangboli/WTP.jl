![WTP LOGO](https://github.com/kangboli/WTP.jl/blob/main/wtp_log.svg?raw=true)
# WTP

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://wtp.kangbo.dev)

## Introduction

WTP.jl provides a **thin layer of abstractions** for processing condense-phase
orbitals. The library aims primarily to save ourselves from the insufferable
agony of indices gymnastics without immolating all the performance.  Our goal
is to make condense-phase codes  more concise and comprehensible with a tolerable
performance sacrifice.

The software design of `WTP.jl` follows the way of a computer scientist (as
opposed to that of a computational scientist). I'm approaching numerical
software through thin layers of composible abstractions instead of a large
`code` with an extensive features and optimization.

Currently, we provide a few functionalities.

1. Nonorthogonal and periodic grids: Support intuitive indexing, iteration, and Fourier transforms. This encapsulates the crystal lattice, the reciprocal lattice, the Brillouin zone, and the homecell.
2. Orbitals: We provide an intuitive interface for efficient inner products, fft, various forms of indexing, folding, and translation.
3. Some IO functionalities with Quantum Espresso's `wfc?.dat` files and
Wannier90's `.amn`/`.mmn` files. WTP reads these files and associates them with a mesh for easy indexing.


