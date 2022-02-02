# WTP

## Introduction

WTP.jl provides a thin layer of abstraction for processing condense-phase
orbitals. The library aims primarily to enable concise and readable code without
sacreficing a significant amount of performance. Our goal is to make
condense-phase codes 10+ times less verbose and 10+ times more scrutable
(subjectively speaking) than their Fortran77 counterparts with less than 20\%
performance loss.

Currently, we provide a few abstractions.

1. Indexable, iterable, and Fourier-transformable non-orthogonal meshes for working with
the crystal lattice, the reciprocal lattice, the Brillouin zone, and the
homecell lattice.
2. Orbitals with an intuitive interface for efficient inner products, fft,
various forms of indexing, translation, and lazy algebras.
3. Replicate Wannier90's evaluation of the center and the spread. A steepest
descent algorithm is included as well (This probably will be moved into a
separate package).
4. Some IO functionalities with Quantum Espresso's `wfc?.dat` files and
Wannier90's `.amn`/`.mmn` files. WTP reads these files and associates them with
a mesh for easy indexing.

To enforce readability, we have a few strict rules specifically against some
practices in scientific code that I'm frightened of 

1. No acronym/abbreviation is allowed for function/type/variable names, unless it is commonly used in published literature. For example, it is acceptable to use "fft" for the fast Fourier transform, but using "bz" for the Brillouin zone, Benzene, or Bezanson is not allowed. Note that single-letter Mathematical symbols are not linguistic acronym/abbreviations, and they can be used. For example, `for k in brillouin_zone` is fine, but `for iknum=ikbgn:ikend; kpnt = bz[:, iknum]` is not allowed.
2. No more than two levels of indentation is allowed for procedural code. This rule does not apply to non-procedural programming. For example, three levels of nested for loops or if statements are not allowed, but nested comprehensions or reductions can be arbitrarily indented. One can substitute many `for/while` loops with more specific types of iteration (broadcasting/mapping/reduction/recursion) and replace `if/then/else` with more specific types of predicates (multiple dispatch/filter/ternary operator/short circuiting). 
3. Do not write comment to explain things. If you think a piece of code requires an explanation, please refactor the code (into functions with appropriate names) until it doesn't or just leave a url to something more illustrative than unicode. Explaining things with comments is mostly useless and sometimes harmful in an expressive language. Also, don't use comment for who did what at what time for what purpose. VCS does that better without polluting the code base. In general, the primary purpose of comments should be sarcarsm.
4. The scope of any variable cannot surpass 30 single lines of code (SLOC), which is less than a page on most monitors with a reasonable font. That is, all functions have to fit in one screen, and global variables are verboten. 