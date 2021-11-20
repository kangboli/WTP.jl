wave_functions_list = wave_functions_from_directory(joinpath(test_3_dir, "benzene.save"))
wannier = wannier_from_save(wave_functions_list);
brillouin_zone = grid(wannier)
wannier_real = ifft(wannier)
u_1 = wannier_real[brillouin_zone[0,0,0]][1]
homecell = grid(u_1)
u_1[homecell[0, 0, 0]]
u_1[homecell[0, 0, 2]]
export_to_gaussian_cube("benzene_1.cube", wannier_real, [1],
[(1, 0., 0., 0., 0.)])
export_to_gaussian_cube("benzene_2.cube", wannier_real, [2],
[(1, 0., 0., 0., 0.)])