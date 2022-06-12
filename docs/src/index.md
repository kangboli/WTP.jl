# WTP

## Introduction

WTP.jl provides a thin layer of low cost abstractions for processing
condense-phase orbitals. The library aims primarily to save ourselves from the
insufferable agony of indices gymnastics without immolating all the performance.

Currently, we provide a few abstractions.

1. Indexable, iterable, and Fourier-transformable non-orthogonal meshes for working with the crystal lattice, the reciprocal lattice, the Brillouin zone, and the homecell lattice.
2. Orbitals with an intuitive interface for efficient inner products, fft, various forms of indexing, translation, and lazy algebras.
3. Some IO functionalities with Quantum Espresso's `wfc?.dat` files and Wannier90's `.amn`/`.mmn` files. WTP reads these files and associates them with a mesh for easy indexing.
4. Implements the convolutional center and spread.

There is a number of taboos for this package.

1. No abbreviations/acronyms. 
2. No more than two levels of indentation.
3. No global variables.
4. No comments on who did what at what time.
5. No functions longer than one page.

