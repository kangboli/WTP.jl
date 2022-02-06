The centers and the spreads of the Wannier functions are significantly more
complicated than the generic mean and the variance in statistics.
Here we try to clarify all the subtleties.

## The ground truth definition

The center of the Wannier function $|\mathbf{0}n\rangle$, which is the $n^{\text{th}}$ wannier function in homecell $\mathbf{0}$,  is defined in the literature as
$$\overline{r}_n = \langle \mathbf{0} n | \mathbf{r} | \mathbf{0} n \rangle}$$.

This is immediately a questionable definition to build upon. 
The definition is a statistical mean, which implicitly assume a hard boundary
conditions. In solid state physics, however, orbitals have periodic boundary
conditions. For orbitals spanning across the boundaries, this definition clearly
makes no sense.

![](cross_boundary.svg)

This mainly causes problems when the number of k-point is small or the orbitals
are diffused. In the case of $\Gamma$-point calculations, this definition breaks
down entirely. 

That being said, Wannier90 works well, but their finite difference formula for
the center approximates this definition only to linear order. In my opinion, their
approximation is significantly more accurate than their definition. 

If we put up with this definition, the spread is

$$\Omega = \sum_n [\langle \mathbf{0}n | \mathbf{r}^2 | \mathbf{0} n \rangle 
 - \langle \mathbf{0} n | \mathbf{r} | \mathbf{0} n \rangle^2 ]$$.

This definition is equally questionable, but again, their is an approximation
that  is more accurate.

## Reciprocal Space Formula

The centers and the spreads are evaluated in the reciprocal space roughly
because many things can be formulated nicely in terms of the nearest neighbor
integrals $M_{mn}^{\mathbf{k}, \mathbf{b}}$.  The reciprocal space formulas are






