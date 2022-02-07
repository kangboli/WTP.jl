# WTP

## Introduction

WTP.jl provides a thin layer of low cost abstractions for processing
condense-phase orbitals. The library aims primarily to save ourselves from the
insufferable agony of indices gymnastics without immolating all the performance.
Our goal is to make condense-phase codes 10+ times more concise and 10+ times
more scrutable (subjectively speaking) than their Fortran77 counterparts with
less than 20\% performance loss.

Currently, we provide a few abstractions.

1. Indexable, iterable, and Fourier-transformable non-orthogonal meshes for working with the crystal lattice, the reciprocal lattice, the Brillouin zone, and the homecell lattice.
2. Orbitals with an intuitive interface for efficient inner products, fft, various forms of indexing, translation, and lazy algebras.
3. Replicate Wannier90's evaluation of the center and the spread. A steepest descent algorithm is included as well (This probably will be moved into a separate package).
4. Some IO functionalities with Quantum Espresso's `wfc?.dat` files and Wannier90's `.amn`/`.mmn` files. WTP reads these files and associates them with a mesh for easy indexing.
5. Implements the convolutional center and spread and teh ILA optimization, but this should probably be moved a separate pacakge along with W90 functionalities.

There is a number of taboos for this package.

1. abbreviations/acronyms. 
2. over two levels of indentation.
3. global variables.
4. comments on who did what at what time.
5. functions longer than one page.
