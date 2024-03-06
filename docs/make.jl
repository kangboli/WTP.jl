using Documenter, WTP, LinearAlgebra

push!(LOAD_PATH,"../src/")
makedocs(sitename="𝑾𝑻𝑷.jl",
modules=[WTP],
checkdocs=:none,
doctest=false,
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
]
)

deploydocs(
    repo = "github.com/kangboli/WTP.jl",
)
