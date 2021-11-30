push!(LOAD_PATH,"../src/")
using Documenter, WTP
makedocs(sitename="𝑾𝑻𝑷.jl",
pages = [
    "Home" => "index.md",
    "Grid" => "grid.md",
    "Grid Vector" => "grid_vector.md",
    "Orbital" => "orbital.md",
    "File IO" => "iopw.md",
    "Wannier" => "wannier.md",
    "Γ-point" => "gamma_point.md",
    "Examples" => "examples.md",
]
)

