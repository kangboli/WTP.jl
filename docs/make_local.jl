using Documenter, WTP, LinearAlgebra

push!(LOAD_PATH,"../src/")
#= DocMeta.setdocmeta!(WTP, :DocTestSetup, :(using WTP; using LinearAlgebra; using UnicodePlots); recursive=true) =#
#
makedocs(sitename="𝑾𝑻𝑷.jl",
modules=[WTP],
checkdocs=:none,
doctest=false,
# strict=true,
#= format=Documenter.HTML(
    assets=["assets/custom.css"],
), =#
pages = [
    "Home" => "index.md",
    "Grid" => "grid.md",
    "Grid Vector" => "grid_vector.md",
    "Conversion of Indices" => "miller_indices.md",
    "Import & Export" => "iopw.md",
    "Function on a Grid" => "on_grid.md",
    "Orbital" => "orbital.md",
    "Orbital Set " => "orbital_set.md",
    # "Center and Spread" => "center_spread.md"
    # "Wannier" => "wannier.md",
    # "Γ-point" => "gamma_point.md",
    # "Examples" => "examples.md",
    # "API References" => "api_references.md",
    # "Internals" => "internals.md"
]
)

