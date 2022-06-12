# WTP

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://wtp.kangbo.dev)
## Introduction

WTP.jl provides a **thin layer of abstractions** for processing condense-phase
orbitals. The library aims primarily to save ourselves from the insufferable
agony of indices gymnastics without immolating all the performance.  Our goal is
to make condense-phase codes  more concise and scrutable with tolerable
performance sacrifice.

The software design of WTP.jl follows the way of a computer scientist (as
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

## Parallelization

I'm not yet sure about how to go about this, but I will almost certainly reject
anything put forth by physicists alone (MPI). One can parallelize code with such
things and get good performance, but the sacrifice is your soul 

- correctness: MPI leaves all the concurrency crap to you.
- flexibility: No more multiple dispatch, functional programming, and meta programming. 
- readability: you are no longer reading one program at a time.
- maintainability: tests now have to be run on a cluster.
- usability: the user have to configure and use MPI to try your package.
- fault tolerance: you have to implement a `restart` mode somewhere.
- modularity: MPI state is global. the whole point of modularity is to avoid sharing states.

I don't think it is the right path forward. 

I always had difficulty talking to system researchers about Fortran+MPI and
opportunities to do something better. The system so absurd to them that they
don't believe it is a real thing, and they consider me as either stupid or
delusional.
