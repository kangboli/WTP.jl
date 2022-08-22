"""
This file contains hacks for interfacing with QE and Wannier90. 
Only code within this file can assume knowledge of the traps
in the code of QE and Wannier90.

The interface has two levels. Level one is a one-to-one translation of
QE/Wannier90 file structure to Julia structs. Level two maps these structs to
programming abstractions in the package. Level one should entirely abstract away
QE/Wannier90, of which level two must be oblivious. This also means 
that the rest of the package, which interface through level two, must be oblivious
of QE/Wannier90.
"""

using FortranFiles

export WFC,
    wave_functions_from_directory,
    orbital_set_from_save,
    load_evc!,
    orbitals_from_wave_functions,
    UNK,
    single_orbital_from_unk,
    wannier_from_unk_dir,
    orbitals_from_unk,
    MMN,
    AMN,
    UMAT,
    wave_function_basis,
    brillouin_zone_from_k_coordinates,
    estimate_sizes,
    i_kpoint_map,
    single_orbital_from_wave_functions,
    dump


function angry_parse(::Type{T}, n) where {T}
    try
        return parse(T, n)
    catch
        return T(0)
    end
end


"""
Parse a single line (string) as an array of primitive types.
"""
parse_line = (line, type) -> (n -> angry_parse(type, n)).(split(strip(line), r"\s+"))

"""
This structure provides an interface to the .wfc files.

The format is documented in: [format of wfc files](https://gitlab.com/QEF/q-e/-/wikis/Developers/Format-of-wfc-files).

The relevant part is copied here for convenience:

- `INTEGER :: ik`  k-point index (1 to number of k-points)
- `REAL(8) :: xk(3)`  k-point coordinates
- `INTEGER :: ispin`  spin index for LSDA case: ispin=1 for spin-up, ispin=2 for spin-down for unpolarized or non-colinear cases, ispin=1 always
- `LOGICAL :: gamma_only`  if .true. write or read only half of the plane waves
- `REAL(8) :: scalef`  scale factor applied to wavefunctions
- `INTEGER :: ngw`  number of plane waves (PW)
- `INTEGER :: igwx`  max number of PW (may be larger than ngw, not sure why)
- `INTEGER :: npol`  number of spin states for PWs: 2 for non-colinear case, 1 otherwise
- `INTEGER :: nbnd`  number of wavefunctions
- `REAL(8) :: b1(3), b2(3), b3(3)`  primitive reciprocal lattice vectors
- `INTEGER :: mill(3,igwx)`  miller indices: `h=mill(1,i), k=mill(2,i), l=mill(3,i)` the i-th PW has wave vector `(k+G)(:)=xk(:)+h*b1(:)+k*b2(:)+ l*b3(:)`
- `COMPLEX(8) :: evc(npol*igwx,nbnd)`  wave functions in the PW basis set The first index runs on PW components, the second index runs on band states.  For non-colinear case, each PW has a spin component first  igwx components have PW with   up spin, second igwx components have PW with down spin

Some variables have been renamed in the package to exorcise the acronyms.

- `ik` -> `i_kpoint`
- `xk` -> `k_coordinates`
- `ispin` -> `i_spin`
- `scalef` -> `scale_factor`
- `ngw` -> `n_planewaves`
- `igwx` -> `max_n_planewaves`
- `npol` -> `n_polerizations`
- `nbnd` -> `n_band`
- `b1, b2, b3` -> `b_1, b_2, b_3`
"""
mutable struct WFC
    filename::String
    i_kpoint::Int64
    k_coordinates::Vector{Float64}
    i_spin::Int64
    gamma_only::Bool
    scale_factor::Float64
    n_planewaves::Int64
    max_n_planewaves::Int64
    n_polerizations::Int64
    n_band::Int64

    b_1::Vector{Float64}
    b_2::Vector{Float64}
    b_3::Vector{Float64}
    miller::Matrix{Int64}
    evc::Union{Matrix{ComplexFxx},Nothing}
end

"""
The basis vectors of the Reciprocal lattice as a matrix.
Each column is the coordinates of a basis vector.
"""
wave_function_basis(w::WFC)::Matrix{Float64} = [w.b_1 w.b_2 w.b_3]

"""
    WFC(wave_function_filename)

Load the metadata of a wave function from a wfc.dat (terrible filenames) file.
This does not load the wave functions.

This is not guaranteed to work when Quantum Espresso is compiled with a Fortran compiler that
is not `gfort`. Fortran does not enforce a standard on the size of it's
primitive types in binary files. 

Example: 

Assume that you are in `docs` directory (change the path otherwise). 

```@jldoctest wfc
julia> path_to_si = "../test/test_data/test_5";
julia> wave_function = WFC(joinpath(path_to_si, "si.save/wfc1.dat"));
julia> wave_function.n_band, wave_function.gamma_only, wave_function.n_planewaves
(4, false, 537)
```
"""
function WFC(wave_function_filename::String)

    data = FortranFile(wave_function_filename)
    i_kpoint, k_coordinates, i_spin, gamma_only, scale_factor =
        read(data, Int32, (Float64, 3), Int32, Int32, Float64)
    gamma_only = Bool(gamma_only)
    n_planewaves, max_n_planewaves, n_polerizations, n_band =
        read(data, Int32, Int32, Int32, Int32)
    b_1, b_2, b_3 = read(data, (Float64, 3), (Float64, 3), (Float64, 3))
    miller = read(data, (Int32, (3, max_n_planewaves)))
    close(data)

    return WFC(
        wave_function_filename,
        i_kpoint,
        k_coordinates,
        i_spin,
        gamma_only,
        scale_factor,
        n_planewaves,
        max_n_planewaves,
        n_polerizations,
        n_band,
        b_1,
        b_2,
        b_3,
        miller,
        nothing,
    )

end

"""
    wave_functions_from_directory(save_dir)

Read all the wave function files in an directory for metadata. This does not
load the wave functions.

The `wfcX.dat` files only gives you the content of a single wave function, but
what you really want is probably reading the entire `.save` directory and
convert them to orbitals on grids. 

```@jldoctest wfc
julia> wave_functions_list = wave_functions_from_directory(joinpath(path_to_si, "si.save"));
julia> length(wave_functions_list)
64
```
"""
function wave_functions_from_directory(save_dir::String)
    is_wave_functions(f::String) =
        startswith(basename(f), "wfc") && endswith(basename(f), ".dat")
    return [WFC(file) for file in readdir(save_dir, join = true) if is_wave_functions(file)]
end

"""
    orbital_set_from_save(
        wave_functions_list,
        [domain_scaling_factor=2,]
    )

Load the wavefunctions and parse them into orbitals.

The `domain_scaling_factor` is possibly either from Nyquist-Shannon sampling theorem.
or from the energy cutoff for the density. 
Settting this to 2 gives the same FFT as QE.

```@jldoctest wfc
julia> u = orbital_set_from_save(wave_functions_list);
julia> length(elements(u))
64
```
"""
function orbital_set_from_save(
    wave_functions_list::AbstractVector{WFC};
    domain_scaling_factor::Integer = 2,
    bands::Vector{Int} = collect(1:wave_functions_list[1].n_band)
)
    k_map, brillouin_zone = i_kpoint_map(wave_functions_list)

    wannier = init_orbital_set(brillouin_zone)
    sizes = Tuple(
        maximum((w) -> estimate_sizes(w, i, domain_scaling_factor), wave_functions_list) for i = 1:3
    )

    for w in wave_functions_list
        k = k_map[w.i_kpoint]
        gauge(wannier)[reset_overflow(k)] = Matrix{Float64}(I, length(bands), length(bands))

        load_evc!(w, bands)
        reciprocal_lattice = make_grid(ReciprocalLattice3D,
            matrix_to_vector3(wave_function_basis(w)),
            size_to_domain(sizes),
        )
        wannier[reset_overflow(k)] = orbitals_from_wave_functions(w, reciprocal_lattice, k, bands)
        w.evc = nothing
    end
    return wannier
end

"""
The returned orbitals will be in the first Brillouin zone.

Create all the orbitals that are in the wave_functions structure.

k can be out of the first Brillouin zone. Images within the first Brillouin zone
will be constructed and returned.
"""
function orbitals_from_wave_functions(
    wave_functions::WFC,
    reciprocal_lattice::ReciprocalLattice,
    k::KPoint,
    bands::Vector{Int}
)
    return [
        single_orbital_from_wave_functions(
            wave_functions,
            reciprocal_lattice,
            k,
            index_band,
            band_order
        ) for (band_order, index_band) in enumerate(bands)
    ]
end


"""
   single_orbital_from_wave_functions(w, l, k, n) 

Use this if you want to load a single orbital. Mostly likely you will want
`orbital_set_from_save` to load them all at once.

Create an orbital on the nth band from a WFC object w. The orbital will be
defined on a reciprocal lattice l.

k is the kpoint, which can be out of the first Brillouin zone.  If k is outside
the first Brillouin zone, we will fold the kpoint into the first Brillouin zone
and translate the wave function in on the reciprocal lattice (phase shift in real lattice).

The orbital is represented as values on the reciprocal lattice.
"""
function single_orbital_from_wave_functions(
    wave_function::WFC,
    reciprocal_lattice::ReciprocalLattice,
    k::KPoint,
    band_index::Int,
    band_order::Int
)
    wave_function.evc !== nothing || error("The wave functions are not loaded.")
    empty_elements = zeros(ComplexFxx, size(reciprocal_lattice)...)
    orbital = UnkBasisOrbital(reciprocal_lattice, empty_elements, reset_overflow(k), band_index)

    i_kpoint!(orbital, wave_function.i_kpoint)

    Threads.@threads for i = 1:wave_function.max_n_planewaves
        coefficients = wave_function.miller[:, i] + overflow(k)
        g = make_grid_vector(reciprocal_lattice, coefficients)
        orbital[g] = wave_function.evc[i, band_order]
    end

    """
    This little function was obnoxiously difficult to 
    figure out. The details are documented in the README.
    """
    function complete_negative_half_for_gamma_trick()
        for g in (m -> reciprocal_lattice[m...]).(eachcol(wave_function.miller))
            orbital[-g] = conj(orbital[g])
        end
    end

    wave_function.gamma_only && complete_negative_half_for_gamma_trick()
    return orbital
end

"""
    load_evc!(wfc, [bands=1:wfc.n_band])

Load the evc (frequency space wave-functions) into wave_functions.
You probably shouldn't have to touch this.
"""
function load_evc!(w::WFC, bands=1:w.n_band)
    data = FortranFile(w.filename)

    # TODO: Demystify this piece of code.
    for _ = 1:4
        read(data)
    end
    w.evc = zeros(ComplexFxx, w.n_polerizations * w.max_n_planewaves, length(bands))
    c = 0
    for i = 1:w.n_band
        band_data = read(data, (ComplexF64, w.n_polerizations * w.max_n_planewaves))
        i ∉ bands && continue
        c += 1
        w.evc[:, c] = band_data
    end
    close(data)
end


"""
    i_kpoint_map(wave_function_list)

Construct a mapping between the mystical i_kpoints and kpoints from a list of
WFC objects. Also figure out the Brillouin zone.  For technical reasons, these
two things are difficult to do separately.

This mapping may be different from what is used for .mmn and .amn files.

```@jldoctest wfc
julia> k_map, brillouin_zone = i_kpoint_map(wave_functions_list);
julia> k_map[1]
GridVector{BrillouinZone3D}:
    coefficients: [0, 0, 0]
```
"""
function i_kpoint_map(wave_functions_list::AbstractVector{WFC})
    k_map = Dict{Int64,KPoint}()

    c = (w -> w.k_coordinates).(wave_functions_list)
    brillouin_zone =
        brillouin_zone_from_k_coordinates(c, wave_function_basis(wave_functions_list[1]))
    for w in wave_functions_list
        k_map[w.i_kpoint] = snap(brillouin_zone, w.k_coordinates)
    end
    return k_map, brillouin_zone
end

"""
Estimate the size of the reciprocal lattice/homecell lattice.

It has to be large enough to hold all the planewaves and some slack space for
phase shifting. A slack of 1 is necessary for correctly representing the wave functions
within the first Brillouin zone because QE/Wannier90 do not use the first Brillouin zone.

Additional slack my be needed for representing wave functions outside the first
Brillouin zone, which occur in finite difference schemes.  The slack factor
depends on the scheme. Gor the finite different scheme in MLWF, a slack of 1
should suffice.
"""
function estimate_sizes(
    wave_functions::WFC,
    i::Int,
    domain_scaling_factor::Number,
    phase_shift_slack::Integer = 1,
)
    domain_scaling_factor == 2 &&
        return 2 * domain_scaling_factor * (max(abs.(wave_functions.miller[i, :])...))
    domain_scaling_factor == 1 ||
        error("Only a domain scaling factor of 1 and 2 are currently supported")
    half_domain = max(abs.(wave_functions.miller[i, :])...)
    return 2 * (half_domain + phase_shift_slack)
end

"""
Construct a BrillouinZone (grid) from a set of kpoints and the reciprocal basis vectors.

This is done by first writing the kpoints in terms of the reciprocal basis.
Then we take the three smallest kpoints besides the gamma point. 
"""
function brillouin_zone_from_k_coordinates(
    k_coordinates::Vector{Vector{Float64}},
    reciprocal_basis::Matrix{Float64},
)
    gamma_point = k_coordinates[argmin(norm.(k_coordinates))]
    offset_kpoints = [k_coordinates[i] - gamma_point for i = 1:length(k_coordinates)]
    in_reciprocal_basis = [reciprocal_basis \ k for k in offset_kpoints]

    """
    Find the size of a 1D grid given the set of kpoint coordinate along that
    direction in the reciprocal basis.
    """
    function find_size(v::Vector{Float64})
        non_zeros = collect(filter((x) -> !isapprox(x, 0, atol = 1e-7), v))
        return length(non_zeros) == 0 ? 1 : round(1 / min(non_zeros...))
    end

    brillouin_sizes = Tuple(find_size((k -> k[i]).(in_reciprocal_basis)) for i = 1:3)
    make_grid(BrillouinZone3D,
        matrix_to_vector3(reciprocal_basis * inv(diagm([brillouin_sizes...]))),
        size_to_domain(brillouin_sizes),
    )
end


"""
The `UNKXXXXX.1` files are produced by `pw2wannier90.x` to show orbitals in real space.

This structure provides an interface to the UNK files.

- `nx, ny, nz`: Number of grid points in x, y, z direction.
- `k`: The kpoint at which the unk orbital is defined.
- `n_band`: The number of bands stored.
- `psi_r`: Each row of this matrix corresponds to a band. Each column corresponds
    to a grid point (not sure which one yet).
- `filename`: The name of the file from which the UNK object is constructed.
"""
struct UNK
    nx::Int64
    ny::Int64
    nz::Int64
    k::Int64
    n_band::Int64
    psi_r::Matrix{ComplexFxx}
    filename::String
end



"""
    UNK(unk_filename)

Load a unk file named unk_filename as a UNK object.

The wave functions are read such that each row is an orbital
to avoid transposing in memory for SCDM.

This takes 10s to load a 1.4G unk file using 12 threads and an NVMe drive.
Since the unk files are too clumsy to cinlude in VCS, we will not provide the files.

Example:

```julia
wave_function = UNK(joinpath(path_to_si, "unk/UNK00001.1"))
```
"""
function UNK(unk_filename::String)
    file = open(unk_filename)
    # TODO: provide an option for using string splitting (slower).
    parse_complex(number) =
        parse(Float32, number[1:20]) + parse(Float32, number[21:40]) * 1im
    n_x, n_y, n_z, k, n_band = parse_line(readline(file), Int64)
    empty_elements = zeros(ComplexFxx, n_band, n_x * n_y * n_z)
    unk = UNK(n_x, n_y, n_z, k, n_band, empty_elements, unk_filename)
    ng = n_x * n_y * n_z
    lines = readlines(file)
    for n = 1:n_band
        Threads.@threads for i = 1:ng
            unk.psi_r[n, i] = parse_complex(lines[ng*(n-1)+i])
        end
    end
    close(file)
    return unk
end

"""
    single_orbital_from_unk(unk, h, k, n)

Use this if you want to load a single orbital. Mostly likely you will want
`wannier_from_unk_dir` to load them all at once.

Create an orbital on the nth band from a UNK structure unk. The orbital will be defined
on a homecell h.

k is the kpoint, which can be out of the first Brillouin zone. 

*The kpoint is not folded and the wave function is not phase shifted*. 
Currently, this is done at a later stage after the fft. 
TODO: Phase shift the orbital.

The orbital is represented as values on the homecell.

"""
function single_orbital_from_unk(unk::UNK, h::HomeCell, k::KPoint, n::Int)
    # TODO: Figure out the layout.
    if (unk.nx != unk.ny) || (unk.nx != unk.nz)
        @warn "UNK files with different nx ny nz dimension has not been
                implemented because pw2wannier90.x has no documentation"
    end

    elements = reshape(unk.psi_r[n, :], (unk.nx, unk.ny, unk.nz))
    return UnkBasisOrbital(h, elements, k, n) |> square_normalize!
end

"""
The returned orbitals will be in the first Brillouin zone.

Create all the orbitals in the unk structure.

A Fourier transform is applied so that orbitals are represented on the
reciprocal lattice.
"""
function orbitals_from_unk(unk::UNK, homecell::HomeCell, k::KPoint)

    """
    Convert the orbital from its image in Wannier brillouin zone to standard brillouin zone.
    """
    function standard_brillouin_zone(orbital, k::KPoint)
        has_overflow(k) || return orbital
        orbital = orbital >> overflow(k)
        kpoint!(orbital, reset_overflow(k))
        return orbital
    end

    function reciprocal_orbital(b)
        o = fft(single_orbital_from_unk(unk, homecell, k, b))
        o = standard_brillouin_zone(o, k)
    end

    reciprocal_orbital.(collect(1:unk.n_band))
end


"""
    wannier_from_unk_dir(
        unk_dir,
        wave_functions_list,
        [domain_scaling_factor=2,]
    )

Load the wavefunctions and parse them into orbitals.
Since the `UNKXXXXX.1` does not carry much information, we need the `wave_functions_list` from 
the `.save` directory as well.

Example: 

```julia
u_real = wannier_from_unk_dir(joinpath(path_to_si, "unk"), wave_functions_list)
```
"""
function wannier_from_unk_dir(
    unk_dir::String,
    wave_functions_list::AbstractVector{WFC},
    domain_scaling_factor::Integer = 2,
)
    k_map, brillouin_zone = i_kpoint_map(wave_functions_list)

    wannier = init_orbital_set(brillouin_zone)
    sizes = Tuple(
        maximum((w) -> estimate_sizes(w, i, domain_scaling_factor), wave_functions_list) for i = 1:3
    )

    for w in wave_functions_list
        k = k_map[w.i_kpoint]
        unk = UNK("$(unk_dir)/UNK$(lpad(w.i_kpoint, 5, "0")).1")
        gauge(wannier)[reset_overflow(k)] = Matrix{Float64}(I, unk.n_band, unk.n_band)

        reciprocal_lattice = make_grid(ReciprocalLattice3D,
            matrix_to_vector3(wave_function_basis(w)),
            size_to_domain(sizes),
        )
        homecell = transform_grid(reciprocal_lattice)
        reciprocal_orbitals = orbitals_from_unk(unk, homecell, k)
        (o -> i_kpoint!(o, w.i_kpoint)).(reciprocal_orbitals)
        wannier[reset_overflow(k)] = reciprocal_orbitals
    end

    return wannier
end

"""
This structure provides an interface to the `MMN` files.

- The first line should tell a joke. The parser won't parse unless tickled.

The format of the `MMN` files can be found in the user guide of `wannier90.x`.

Γ-point:

`wannier90.x` uses MMN files also for the X, Y, Z matrices of gamma point
calculations.  This is nowhere documented.

The ``g``-vectors (taken from wannier90 user guide):

The last three integers specify the ``G`` vector, in reciprocal lattice units,
that brings the k-point specified by the second integer, and that thus lives inside
the first Brillouin zone, to the actual ``k + b`` we need.
"""
struct MMN
    n_band::Int64
    n_kpoint::Int64
    n_neighbor::Int64
    integrals::Dict{Pair{Integer,Integer},Matrix{ComplexFxx}}
    translations::Dict{Pair{Integer,Integer},Tuple{Integer,Integer,Integer}}
end

"""
    MMN(mmn_filename, [ad_hoc_gamma = false])

Load a MMN file. Set `ad_hoc_gamma = true` if it is a gamma-only calculation
so that the gamma-trick shenanigans are dealt with.


```@jldoctest wfc
julia> mmn = MMN(joinpath(path_to_si, "output/pw2wan/si.mmn"))
julia> mmn.n_kpoint, mmn.n_neighbor
(64, 8)
```
"""
function MMN(mmn_filename::String, ad_hoc_gamma = false)
    file = open(mmn_filename)
    parse_complex(number) =
        parse(Float64, number[1:20]) + parse(Float64, number[21:end]) * 1im
    _the_first_line_is_comment = readline(file)
    n_band, n_kpoint, n_neighbor = parse_line(readline(file), Int64)
    mmn = MMN(n_band, n_kpoint, n_neighbor, Dict(), Dict())
    lines = readlines(file)

    function construct_matrix_from(start::Integer)
        mat = zeros(ComplexF64, n_band, n_band)
        Threads.@threads for j = 1:n_band 
            mat[j, :] = (i->parse_complex(lines[start+(i-1)*n_band+j])).(1:n_band)
        end
        return mat
    end

    r_1, r_2 = n_neighbor * (n_band^2 + 1), (n_band^2 + 1)
    for l_1 = 1:n_kpoint
        for l_2 = 1:n_neighbor
            start = (l_1 - 1) * r_1 + (l_2 - 1) * r_2 + 1
            i_kpoint, i_neighbor, g_x, g_y, g_z = parse_line(lines[start], Int64)
            ad_hoc_gamma && (i_neighbor = l_2)
            mmn.translations[i_kpoint=>i_neighbor] = (g_x, g_y, g_z)
            mmn.integrals[i_kpoint=>i_neighbor] = construct_matrix_from(start)
        end
    end
    close(file)
    return mmn
end


"""
The `seedname_u.mat` file.

From Wannier90 User guide.

OUTPUT. Written if write_u_matrices = .TRUE.. The first line gives the date and
time at which the file was created. The second line states the number of
kpoints num_kpts and the number of wannier functions num_wann twice. The third
line is empty.

"""
struct UMAT
    n_band::Int64
    n_kpoint::Int64
    u_matrices::Dict{Vector{Float64}, Matrix{ComplexFxx}}
end


"""
    UMAT(u_mat_filename)

Create a UMAT object from a file
"""
function UMAT(u_mat_filename::String)
    file = open(u_mat_filename)
    parse_complex(number) =
        parse(Float64, number[1:16]) + parse(Float64, number[17:end]) * 1im
    _the_first_line_is_comment = readline(file)
    n_kpoint, n_band, _ = parse_line(readline(file), Int64)
    umat = UMAT(n_band, n_kpoint, Dict())
    lines = readlines(file)

    function construct_matrix_from(start::Integer)
        mat = zeros(ComplexF64, n_band, n_band)
        Threads.@threads for j = 1:n_band 
            mat[j, :] = (i->parse_complex(lines[start+(i-1)*n_band+j])).(1:n_band)
        end
        return mat
    end

    r = n_band^2 + 2
    for l = 1:n_kpoint
        start = (l-1) * r + 1
        k_fraction = parse_line(lines[start+1], Float64)
        umat.u_matrices[k_fraction] = construct_matrix_from(start+1)
    end
    return umat
end


"""
The `.amn` files are used by Wannier90 for storing the gauge transform.
"""
struct AMN
    n_band::Int64
    n_kpoint::Int64
    n_wannier::Int64
    gauge::Dict{Int64,Matrix{ComplexFxx}}
end

"""
    AMN(amn_filename)

Create an AMN object from a file.

```@jldoctest wfc
julia> amn = AMN(joinpath(path_to_si, "output/pw2wan/si.amn"));
julia> amn.n_kpoint
64
```
"""
function AMN(amn_filename::String)
    file = open(amn_filename)
    _ = readline(file) # Skip the date.
    n_band, n_kpoint, n_wannier = parse_line(readline(file), Int64)
    amn = AMN(n_band, n_kpoint, n_wannier, Dict())
    lines = readlines(file)
    r_1 = n_band^2
    c = ReentrantLock()

    function update_gauge(i_kpoint, m, n, value)
        lock(c) do
            haskey(amn.gauge, i_kpoint) || (amn.gauge[i_kpoint] = zeros(ComplexFxx, (n_band, n_wannier)))
            amn.gauge[i_kpoint][m, n] = value
        end
    end

    for l_1 = 1:n_kpoint
        Threads.@threads for l_2 = 1:r_1
            line_number = (l_1 - 1) * r_1 + l_2
            m, n, i_kpoint, real_part, complex_part =
                split(strip(lines[line_number]), r"\s+")
            m, n, i_kpoint = (s -> parse(Int64, s)).([m, n, i_kpoint])
            real_part, complex_part = (s -> parse(Float64, s)).([real_part, complex_part])
            update_gauge(i_kpoint, m, n, real_part + complex_part * 1im) 
        end
    end
    close(file)
    return amn
end

"""
    dump(amn, amn_filename, scdm_parameters=(0, 0))

Write an `AMN` object to disk as a `.amn` file.
"""
function dump(amn::AMN, amn_filename::String; scdm_parameters=(0, 0))
    file = open(amn_filename, "w+")
    write(file, "( ͡° ͜ʖ ͡°)\n")
    if scdm_parameters != (0, 0)
        
        @printf file "%8i%8i%8i  %10.8f%10.8f\n" amn.n_band amn.n_kpoint amn.n_wannier scdm_parameters...
    else
        
        @printf file "%12i%12i%12i\n" amn.n_band amn.n_kpoint amn.n_wannier
    end

    for i = 1:amn.n_kpoint
        for (r, c) in Iterators.product(1:amn.n_band, 1:amn.n_wannier)
            value = amn.gauge[i][r, c]
            @printf file "%5i%5i%5i%18.12f%18.12f\n" r c i real(value) imag(value)
        end
    end

    close(file)
end

"""
    Gauge(grid, amn, k_map, [orthonormalization = true])

Create a gauge from an `AMN` object and a `k_map`.

An `AMN` object also just stores the raw data, and we need to convert this to a
`Gauge` to use it. Like `.mmn` files, we need a mapping of k-points as well as
the Brillouin zone to construct the gauge, which can then be indexed with
k-points.

```@jldoctest wfc
julia> U = Gauge(brillouin_zone, amn, k_map);
julia> U[brillouin_zone[0, 0, 0]][:, :]
4×4 Matrix{ComplexF64}:
  -0.424021+0.0687693im  -0.530967+0.0861146im  -0.469845+0.0762009im  -0.540275+0.0876234im
 0.00775601-0.282049im     -0.3875-0.232102im    0.714604+0.092279im   -0.246713+0.369214im
   0.837758+0.116721im   -0.437746+0.134205im   -0.142503-0.20453im    -0.103362-0.045631im
  0.0694477+0.124822im    0.516605+0.173543im    0.152318-0.411003im   -0.694672+0.0889086im
```
"""
function Gauge(grid::Grid, amn::AMN, k_map::Dict{Int,KPoint}, orthonormalization = true)
    gauge = Gauge(grid)

    orthonormalize(A::AbstractMatrix) =
        let (U, _, V) = svd(A)
            U * adjoint(V)
        end

    for (i_kpoint, u) in amn.gauge
        gauge[k_map[i_kpoint]] = orthonormalization ? orthonormalize(u) : u
    end

    return gauge
end

"""
    Gauge(grid, umat)

Create a gauge from a `UMAT` object and a `grid` (brillouin_zone).
"""
function Gauge(grid::Grid, umat::UMAT)
    gauge = Gauge(grid)

    find_second_smallest(x) = sort(unique(x))[2]
    k_unit_fraction = abs.(hcat(keys(umat.u_matrices)...)) |> eachrow  .|> find_second_smallest

    for (k, u) in umat.u_matrices
        gauge[grid[Int.(k ./ k_unit_fraction)...]] = u
    end
    return gauge
end


"""
    AMN(U, k_map)

Construct an `AMN` object from a gauge `U`. 
"""
function AMN(U::Gauge, k_map::Dict{Int, KPoint})
    n_band, n_wannier = size(elements(U)[1])
    
    brillouin_zone = grid(U)
    amn = AMN(n_band, length(brillouin_zone), n_wannier, Dict{Int, Matrix{ComplexFxx}}())
    for (i_kpoint, k) in k_map
        amn.gauge[i_kpoint] = U[k]
    end

    return amn
end
