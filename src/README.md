# WTP

WTP.jl (not an acronym; stands for nothing) is a basic and user-friendly library
for processing Condense-Phase orbitals. The goal is primarily to promote
readable code in our localization package by providing a thin layer of
abstractions. The hope is that we won't have to twist our mortal parts in the
unbearable agony of indices gymnastics. Currently, we provide

1. Indexable, iterable, and Fourier-transformable non-orthogonal meshes for working with
the crystal lattice, the reciprocal lattice, the Brillouin zone, and the
homecell lattice.
2. Orbitals with an intuitive interface for efficient inner products, fft,
various forms of indexing, translation, and lazy algebras.
3. Wannier90's finite difference scheme for evaluation the center and the spread
of Wannier orbitals. `gamma_only` is also supported. An equality efficient higher
order scheme is in development.
4. Some IO functionalities with Quantum Espresso's `wfc?.dat` files and
Wannier90's `.amn`/`.mmn` files. WTP reads these files and associates them with
a mesh for easy indexing.

Some features are considered.

1. A cache for integrals. This can lead to more transparent code, where formulas
are written in terms of orbitals instead of obscure entries of various
multi-indexed integral tables.
