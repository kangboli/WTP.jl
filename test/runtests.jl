using WTP
using Test
using LinearAlgebra
using Profile

const test_1_dir = "./test/test_data/test_1"
const test_2_dir = "./test/test_data/test_2"
const test_3_dir = "./test/test_data/test_3"

include("test_grid.jl")
include("test_grid_vector.jl")
include("test_fft.jl")
include("test_wfc.jl")
include("test_wannier_load.jl")
include("test_orbital.jl")

include("test_linear_combination.jl")
include("test_integral_table.jl"))