using SnoopCompile

precompile(wave_functions_from_directory, (String, ))
precompile(orbital_set_from_save, (Vector{WFC}, Integer, Vector{Int}, ))



