# WTP

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://wtp.kangbo.dev)

username: matrix

passcode: 6210

## Introduction

WTP.jl provides a **thin layer of abstractions** for processing condense-phase
orbitals. The library aims primarily to save ourselves from the insufferable
agony of indices gymnastics without immolating all the performance.  Our goal
is to make condense-phase codes  more concise and comprehensible with a tolerable
performance sacrifice.

The software design of `WTP.jl` follows the way of a computer scientist (as
opposed to that of a computational scientist). I'm approaching numerical
software through thin layers of composible abstractions instead of a piece of
`code` with an extensive number of features and massive parallelization.

Currently, we provide a few functionalities.

1. Indexable, iterable, and Fourier-transformable non-orthogonal meshes for working with
the crystal lattice, the reciprocal lattice, the Brillouin zone, and the
homecell lattice.
2. Orbitals with an intuitive interface for efficient inner products, fft,
various forms of indexing, translation, and lazy algebras.
3. Some IO functionalities with Quantum Espresso's `wfc?.dat` files and
Wannier90's `.amn`/`.mmn` files. WTP reads these files and associates them with
a mesh for easy indexing.
4. Implements the convolutional center and spread.

