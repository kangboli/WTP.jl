# WTP

## Introduction

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

This package prioritizes code readability and documentation coverage.
Performance is a secondary concern, so we will not compromise the clarity of the
code in exchange for a marginal to moderate performance gain, which is even more negligible when placed in any Quantum Chemistry pipeline. 

We will also not provide more features than a single person can maintain.  We
intend to keep things small and isolated instead of providing  a large set of
features and involving a community. Relevant features such as plotting and
interpolating should have their own separate packages that use this package.

To enforce readability, we have a few strict rules specifically against a few
questionable practices in scientific code

1. No acronym/abbreviation is allowed for function/type/variable names, unless it is commonly used in published literature. For example, it is acceptable to use "fft" for the fast Fourier transform, but using "bz" for the Brillouin zone, Benzene, or Bezanson is not allowed. Note that single-letter Mathematical symbols are not linguistic acronym/abbreviations, and they can be used. For example, `for k in brillouin_zone` is fine, but `for iknum=ikbgn:ikend; kpnt = bz[:, iknum]` is not allowed.
3. No more than two levels of indentation is allowed for procedural code. This rule does not apply to functional programming or meta programming. For example, three levels of nested for loops or if statements are not allowed, but nested functional expressions such as comprehensions or reductions can be arbitrarily indented. 
4. Comments will be stripped. Do not use comment for documentation or (Version Control System) VCS purposes. If you think a piece of code requires an explanation, refactor the code until it doesn't. If you want to do anything else with comments,  do that thing with git instead.
6. The scope of a variable cannot surpass 30 single lines of code (SLOC), which is less than a page on most monitors with a reasonable font. That is, all functions have to fit in one page, and global variables are verboten. Try to break large functions down to smaller ones. If the list of argument grows too long, make a `struct`.
