## Multi-demensional Fill

Do not use `fill` in this package. As of Julia 1.7, it has a memory bug when
filling a 3D array of objects.

## Test dependencies

Test dependencies are not separated from the package dependences since Julia has
no acceptible solution to test dependencies as of 1.7.
A few solutions that do not work.

1. Use `[extras]` and `[targets]` in `Project.toml`. This works only with `Pkg.test()`. You cannot just run the test code in the repl.
2. Have a `test/Project.toml`. This file will not have the package (WTP) in it.
Our package package is only added upon `Pkg.test()`, so we cannot use WTP outside
`Pkg.test()`.
3. If we manully include WTP in `test/Project.toml`, we would be able to use WTP outside `Pkg.test()`, but `Pkg.test()`  would break with an error message `cannot merge packages`.

What I decided to do is to not use test dependencies at all, which only involve
`Test` and `Profile` for us. The downside is that we have more dependencies than we
actually need, but these two are not too much of an issue since they are part of Julia core.