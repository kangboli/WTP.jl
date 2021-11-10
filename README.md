# WTP

## The Brillouin Zone 

Fortran does not do negative indices. Instead of implementing an abstraction for
negative indices, here is what people came up with

- For the input files, the Brillouin zone is centered at the origin and is defined on $[-k_{max}, k_{max}]$. This is the convention in the literature.
- For the output files, the Brillouin zone starts from the origin and is defined on ($[0, 2 k_{max}]$), not the theoretical convention.
- For the $u_{nk}$ orbitals, $[0, 2k_{max}]$ is stored.
- The indices are $0, 1, \ldots, 2 n-1$.

We will use a conventional and sane Brillouin zone by implementing an abstraction for negative indices.

*The Brillouin zone is defined on $[-k_{max}, k_{max}]$ and index by $[-n_{max}+1, n_{max}]$*. (index starts at $1$ in Julia)


## The reciprocal lattice

The reciprocal lattice is $[-nG, nG]$. The lattice points are indexed
by $-n, \ldots, n-1$. For us, the lattice is $[-nG, nG]$, but the lattice
points are indexed by $-n+1, \ldots, n$.

## FFT and Normalization

We use the same convention as QE in terms of FFT. Forward FFT takes 
an orbital defined on the Homecell and produces an orbital on the reciprocal lattice.
The inverse FFT takes an orbital defined on the Reciprocal lattice and produces 
and orbital on the Homecell.

```
                          FFT         
                       ----------->   
Orbitals on Homecell      iFFT        Orbitals on the Reciprocal Lattice
                       <-----------                                     

```

The normalization convention of QE seems to be 

1. Orbitals on the reciprocal lattice normalize to the number of grid points $N_g$.
2. Orbitals on the homecell normalize to $\sqrt{N_g}$.
3. The norm is preserved after applying an inverse Fourier transform and a forward Fourier transfrom `ifft(fft(orbital))`.

I'm not sure if there is reason for this convention, but I find it difficult, under this convention, to 

1. keep operations such as inner products efficient. It is cheaper to keep the
wave functions normalized than adding a normalization step to inner products.
2. keep track of the normalization factor. Many stupid mistakes can be avoided if things normalize to $1$.

For this program, we will normalize the underlying arrays of orbitals to $1$.


## The MMN File

Fortran is a langauge that predates grammar. Instead of learning some grammar, here is what people came up with.

### The Header

The first line is the date of creation (cringe). 
The second line is a list of three integers. 

```
n_band      n_kpoint       n_neighbor
```

### Neighbors

Each k-point has `n_neighbor` neighbors. For each pair of k-points and neighbors,
there is a line of five integers

```
i_kpoint    i_neighbor      gx      gy      gz
```

Someone had the idea of representing every k-point as an integer without 
providing a function for recovering the k-point. This is not a new idea. In computer security, this is
called an **encryption**. The mapping is stored in the `.wave_functions` files, which are in binary format.

TODO: `gx gy gz` 

### The Matrices 

The matices are stored in column major order.
So line `(i-1) * n_band + j` stores row `j` column `i`, which
is the inner product of band `j` of `k` and band `i` of `k+b`.
In other words, $\langle u_{j k} | u_{i k+b} \rangle$ is at `(i-1) * n_band + j`.

## Indexing 3D Objects.

**3D objects should be indexed by a 3D vector**.

Conversion of 3D indices and 1D indices should be abstracted away. There should
never be a varaible like `ik` in QE. You can't blame this on Fortran since Fortran
has first class support for multi-dimensional arrays from the start.


## Programming Guide and Rants

Take a look at the programming practices of Quantum Espresso, and avoid them all.

### Naming

**Abbreviations and acronyms are forbidden**.  

Better be verbose than cryptic. If a name is too long, introduce a **single**
letter alias within a small scope (20 SLoC) so that the full name is visible on the same page.
Never abbreviate anything with a hideous sequence of consonants. It might make you feel
like a hacker, but someone with dyslexia wonders what the heck.

The only special cases are 

1. The acronym is common knowledge. For example, "fft" (fast Fourier transform) is not a problem, but "bz" (Brillouin Zone) is not acceptable. It may as well mean Benzene, buzz, or Benz.

when interfacing with QE/Wannier90. We reuse the
mystical variable names to make the correspondence explicit.


### Global Variables.

**Mutable global variables are forbidden**. 

All mutable variables should be local. Global variables are allowed only for
constants. It is universally agreed upon that using mutable global variables is
a bad practice.  QE not only uses them, it manages to go to the next level and
share them across different programs!


### Indentation

**No more than two levels of indentation is allowed for control flow**.

Code becomes thoroughly incomprehensible when deeply nested. Indentations in QE 
can wrap around the line if you are in protrait mode. Generally, it is never neccessary to
indent more than 2 levels. To reduce the level of indentation, there is a handful of tricks:

1. Abstract any indented peice of code into a function.
2. Use `&&`, `||`, and `? :` expressions instead of `if` statement for early `returns`/`continue`.
3. Use broadcasting/comprehension for parallel code, and recursions for complex iterative code.
    No `for`/`while` loops shall be nested.

### Documentation

**Comments will be ignored**.

Avoid using comments in the source code for documentation purposes. If a piece
of code requires documentation, refactor it into a function and document it with
an appropriate function name and a doc string with usage. Using comments is
otherwise fine in scripts and tests.

