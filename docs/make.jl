using Documenter, WTP, LinearAlgebra

push!(LOAD_PATH,"../src/")
DocMeta.setdocmeta!(WTP, :DocTestSetup, :(using WTP; using LinearAlgebra; using UnicodePlots); recursive=true)
makedocs(sitename="𝑾𝑻𝑷.jl",
modules=[WTP],
checkdocs=:none,
pages = [
    "Home" => "index.md",
    "Grid" => "grid.md",
    "Grid Vector" => "grid_vector.md",
    "Indices" => "miller_indices.md",
    "Import/Export" => "iopw.md",
    "OnGrid" => "on_grid.md",
    "Orbital" => "orbital.md",
    "Orbital Set " => "orbital_set.md"
    # "Wannier" => "wannier.md",
    # "Γ-point" => "gamma_point.md",
    # "Examples" => "examples.md",
    # "API References" => "api_references.md",
    # "Internals" => "internals.md"
]
)

