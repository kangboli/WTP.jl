# WTP Design Documentation

## Conventions for Grids

### Standard Form

A grid is in the standard form if

+ The grid is odd and the indices are $-N, \ldots,  N$. This type grids are currently not tested.
+ The grid is even and the indices are $-N, \ldots, N-1$. We will denote this grid $\mathbb{G}_N$.

This convention is not natural to Julia, in which an index starts at $1$, but it seems 
to be more consistent with existing condense phase code, although there doesn't appear to 
be a unified convention for grids in that literature.

### Grid Center

The center of a grid that is in the standard form is $0$. If a grid is not in the standard form, 
it can be obtained by applying an unique translation to some grid in standard form. The translation
is defined as the center. 

For example, an even grid with the indices $-N+1, \ldots, N$ can be obtained by translating
$\mathbb{G}_N$ to the right by $1$. So center of the grid is $1$.

### Array Storage

Any object defined on the grid is stored as an (most likely 3D) array. The
layout of the array, when the grid is in the standard form, is as follows

- An non-negative grid index $i$ corresponds to the array index $i+1$. For example, the center
of the grid corresponds to the first slot of the array.
- A negative grid index $i$ corresponds to $i + s + 1$, where $s$ is the size of the grid. i.e. Negative grid indices wrap around the array.

There are more complications when the grid is not in the standard form. Non-standard grid is 
not yet tested or documented.

On-grid objects are stored as such primarily because this is what `fft` wants.

### Mappings of Indices

An on-grid object can be indexed in a few different ways depending on the usage

1. Indexed by a grid vector. This is
notationally the most convenient way to index an on-grid object such as a
function defined on the reciprocal lattice.
2. Indexed by three array indices. This is useful for manipulating the underlying array
for low-level operations. Note that these indices are very different from the coefficients
of a grid vector.
3. Indexed by a single number. This is useful if one needs to treat an on-grid object 
as a 1D vector to apply standard matrix/vector algorithms (such as QRCP).

In principle, all the exported functions should only accept and return grid
vectors except for the functions that convert the indices. This is because grid
vectors can be converted to the other forms without a context, whereas the
opposite is not true.  Array indices and single index can be used strictly
inside top-level functions, where the context is available.

Some utility functions are provided for easy conversion amongst the three.

- `three_to_one` converts a single index to array indices. 
- `one_two_three` converts array indices to a single index.
- `miller_to_standard` converts a grid vector (or just its coefficient) to array indices. 
- `standard_to_miller` converts a set of array indices to a grid vector.

## The Brillouin Zone 

Fortran does not do negative indices. Instead of an abstraction for
negative indices, here is what we got in some Fortran code

- For the input files, the Brillouin zone is centered at the origin $[-k_{max}, k_{max}]$. This is the convention in the literature.
- For the output files and internal use, the Brillouin zone starts from the origin and is defined on $[0, 2 k_{max}]$. This is irritating.
- The $u_{nk}$ orbitals stored correspond to $[0, 2k_{max}]$.
- The indices are $0, 1, \ldots, 2 n-1$.

We use our unified abstraction of grid indices, which matches the convention in the literature.
To reiterate,

*The Brillouin zone is defined on $[-k_{max}, k_{max}]$ and index by $[-n_{max}, n_{max}-1]$*. 

This difference in convention turns out to be quite irritating.  For orbitals
outside the first Brillouin zone, QE store $u_{nk}$ with $k$ outside the first
Brillouin zone. To get its image within the first Brillouin zone.  We need to
phase shift it. This means a translation (circular shift) in the frequency
space. However, if our grid is just large enough to contain the orbital, a
circular shift would give a wrong result. Because of this, we have to zero-pad
our orbitals in the reciprocal lattice and introduce another difference from
QE.


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
3. The norm is preserved after applying an inverse Fourier transform and a forward Fourier transform `ifft(fft(orbital))`.

I'm not sure if there is reason for this convention, but I find it difficult, under this convention, to 

1. keep operations such as inner products efficient. It is cheaper to keep the
wave functions normalized than adding a normalization step to inner products.
2. keep track of the normalization factor. Many stupid mistakes can be avoided if things normalize to $1$.

For this program, we will normalize the underlying arrays of orbitals to $1$.

When we Fourier transform an orbital, we transform not only the wave function,
but also the grid on which it is defined. A homecell will be transformed into a
reciprocal lattice; a Brillouin zone will be transformed into a Crystal lattice;
and vice versa.

## File IO

### The MMN File

Note that the following description does not apply when the gamma trick is used.
See the gamma trick section for additional information.

#### The Header

The first line is the date of creation (cringe). 
The second line is a list of three integers. 

```
n_band      n_kpoint       n_neighbor
```

#### Neighbors

Each k-point has `n_neighbor` neighbors. For each pair of k-points and neighbors,
there is a line of five integers

```
i_kpoint    i_neighbor      gx      gy      gz
```

Someone had the idea of representing every k-point as an integer without 
providing a function for recovering the k-point. This is not a new idea. In computer security, this is
called an **encryption**. The mapping is stored in the `.wave_functions` files, which are in binary format.

TODO: `gx gy gz` 

#### The Matrices 

The matices are stored in column major order.
So line `(i-1) * n_band + j` stores row `j` column `i`, which
is the inner product of band `j` of `k` and band `i` of `k+b`.
In other words, $\langle u_{j k} | u_{i k+b} \rangle$ is at `(i-1) * n_band + j`.

### The Gamma Trick

TODO: Figure out why there is only factor of 2 instead of 8 reduction.

The gamma trick is scarcely documented and confused me very much.
In my opinion, it really shouldn't be a stand-along trick, which 
pollutes the code base with an exception case all over the place.
It should be abstracted away as part of the storage layer.

#### The Trick.

The trick is that only half (why not 1/8) of the wave-function needs to be 
stored when the calculation is `gamma_only`. The question that is nowhere
documented in the code base is what "half" means. After much reverse engineering,
it turns out that their (positive) half can be characterized as

```
i * n_y * n_z + j * n_z + k >= 0
```

Here `i`, `j`, `k` are the array indices and `n_y`, `n_z` are the array size.
This design is painful to work with since people have to be constantly 
aware of the trick and make a special case every step of the way.
Thus, we will just implement this in one place as a storage optimization that 
is abstracted away from the programming interface.

## Programming Guide and Rants

Take a look at what the Physicists do, and you know what to avoid.

### Naming

> Abbreviations and acronyms are forbidden.

Better be verbose than cryptic. If a name is too long, introduce a **single**
letter alias within a small scope (20 SLoC) so that the full name is visible on
the same page. Especially don't abbreviate anything with a hideous sequence of
consonants.

The only special cases are 

1. The acronym is common knowledge. For example, "fft" (fast Fourier transform) is not a problem, but "bz" (Brillouin Zone) is not acceptable. It may very well mean Benzene, buzz, booz, bizarre, or Benz.
2. When interfacing with other packages. One can reuse the mystical variable names to make the correspondence explicit.

Also, don't spell the same word two ways, especially when they are next to each other. 
Write "localization" instead of "localisation".

### Global Variables.

> Mutable global variables are forbidden. 

All mutable variables should be local. Global variables are allowed only for
constants. If one wants to share a state, pass it as an argument.


### Indentation

> More than two levels of indentation is considered excessive.

Code becomes thoroughly incomprehensible when deeply nested. Indentations in
older electronic structure code can wrap around the monitor if you are in
portrait mode. Theoretically, it is never necessary to indent more than 2 levels.
To reduce the level of indentation, there is a handful of tricks:

1. Abstract any indented piece of code into a function.
2. Use `&&`, `||`, and `? :` expressions instead of `if` statement for early `return`/`continue`.
3. Use broadcasting/comprehension for parallel code, and recursions for complex
   iterative code instead of nested `for`/`while` loops.

### Comments

> If a piece of code is obscure and needs comments, refactor it until it does not.

Avoid using comments in the source code for documentation purposes. If the code
itself cannot be written in a comprehensible way, it cannot be explained in
some figurative allegories that are to be interpreted. Comments also don't
update themselves when you update the code, after which the comments become
lies and go on to gaslight future developers.


### Parallelization

> It won't be MPI.

I'm not yet sure about how to go about this, but I will almost certainly not
use MPI. One can parallelize code with such things and get good performance,
but the sacrifice is everything else

- correctness: MPI leaves all the concurrency crap to you.
- flexibility: No more multiple dispatch, functional programming, and meta programming. 
- readability: you are no longer reading one program at a time.
- maintainability: tests now have to be run on a cluster.
- usability: the user have to configure and use MPI to try your package.
- fault tolerance: you have to implement a `restart` mode somewhere.
- modularity: MPI state is global. the whole point of modularity is to avoid sharing states.

