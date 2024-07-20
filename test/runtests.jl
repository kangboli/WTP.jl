current_dir = pwd()
test_dir = endswith(current_dir, "WTP.jl") ? joinpath(current_dir, "test") : current_dir
test_1_dir = "$(test_dir)/scdm_dataset/test_1"
test_2_dir = "$(test_dir)/scdm_dataset/test_2"
test_3_dir = "$(test_dir)/scdm_dataset/Benzene"
test_6_dir = "$(test_dir)/scdm_dataset/test_6"

include("$(test_dir)/test_grid.jl")
include("$(test_dir)/test_grid_vector.jl")
include("$(test_dir)/test_linear_combination.jl")
include("$(test_dir)/test_miller_indices.jl")
include("$(test_dir)/test_wfc.jl")
include("$(test_dir)/test_io.jl")
include("$(test_dir)/test_wannier_load.jl")
include("$(test_dir)/test_fft.jl")
# include("test_integral_table.jl")
# include("test_gauge_transform.jl")
include("$(test_dir)/test_translation.jl")
# include("test_si_center.jl")
# include("test_gamma.jl")
include("$(test_dir)/test_convolution.jl")
# include("test_w90_gradient.jl")
# include("test_w90_gradient_descent.jl")
# include("test_truncated_convolution.jl")
# include("test_ila_gradient_descent.jl")
# ILA gradient is not tested since there is no direct comparison we can draw.
