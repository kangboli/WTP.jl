"""
Horrible code for dealing with horrible file formats.

This file contains hacks for interfacing with QE and Wannier90. 
Only code within this file can assume knowledge of the shenanigans
in the code of QE and Wannier90.
"""

using FortranFiles

export WFC,
    wave_functions_from_directory,
    wannier_from_save,
    load_evc!,
    orbitals_from_wave_functions,
    UNK,
    single_orbital_from_unk,
    wannier_from_unk_dir,
    orbitals_from_unk,
    MMN,
    AMN,
    wave_function_basis,
    brillouin_zone_from_k_coordinates,
    estimate_sizes,
    i_kpoint_map


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

Load a wave function and its metadata from a wfc.dat (terrible filenames) file.

This may not work when Quantum Espresso is compiled with a different Fortran
compiler.  Fortran does not enforce a standard on the size of it's primitive
types in binary files.  The ways to get consistent behaviors are

1. Stay away from Fortran binary files. 
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
    i_kpoint_map(wave_function_list)

Construct a mapping between the mystical i_kpoints and kpoints
from a list of WFC objects.

This mapping may be different from what is used for .mmn and .amn files.
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
   single_orbital_from_wave_functions(w, l, k, n) 

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
    n::Int,
)
    wave_function.evc !== nothing || error("The wave functions are not loaded.")
    empty_elements = zeros(ComplexFxx, size(reciprocal_lattice)...)
    orbital = UnkBasisOrbital(reciprocal_lattice, empty_elements, reset_overflow(k), n)

    i_kpoint!(orbital, wave_function.i_kpoint)

    Threads.@threads for i = 1:wave_function.max_n_planewaves
        coefficients = wave_function.miller[:, i] + overflow(k)
        g = make_grid_vector(reciprocal_lattice, coefficients)
        orbital[g] = wave_function.evc[i, n]
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
    single_orbital_from_unk(unk, h, k, n)

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
    return UnkBasisOrbital(h, elements, k, n) |> wtp_normalize!
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
)
    return [
        single_orbital_from_wave_functions(
            wave_functions,
            reciprocal_lattice,
            k,
            index_band,
        ) for index_band = 1:wave_functions.n_band
    ]
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
Load the evc (frequency space wave-functions) into wave_functions.
"""
function load_evc!(w::WFC)
    data = FortranFile(w.filename)

    # TODO: Demystify this piece of code.
    for _ = 1:4
        read(data)
    end
    w.evc = zeros(ComplexFxx, w.n_polerizations * w.max_n_planewaves, w.n_band)
    for i = 1:w.n_band
        w.evc[:, i] = read(data, (ComplexF64, w.n_polerizations * w.max_n_planewaves))
    end
    close(data)
end


function wave_functions_from_directory(save_dir::String)
    is_wave_functions(f::String) =
        startswith(basename(f), "wfc") && endswith(basename(f), ".dat")
    return [WFC(file) for file in readdir(save_dir, join = true) if is_wave_functions(file)]
end



"""

The domain_scaling_factor is possibly either from Nyquist-Shannon sampling theorem.
or from the energy cutoff for the density. 
Settting this to 2 gives the same FFT as QE.
"""

function wannier_from_save(
    wave_functions_list::AbstractVector{WFC},
    domain_scaling_factor::Integer = 2,
)
    k_map, brillouin_zone = i_kpoint_map(wave_functions_list)

    wannier = init_wannier(brillouin_zone)
    sizes = Tuple(
        maximum((w) -> estimate_sizes(w, i, domain_scaling_factor), wave_functions_list) for i = 1:3
    )

    for w in wave_functions_list
        k = k_map[w.i_kpoint]
        gauge(wannier)[reset_overflow(k)] = Matrix{Float64}(I, w.n_band, w.n_band)

        load_evc!(w)
        reciprocal_lattice = make_grid(ReciprocalLattice3D,
            matrix_to_vector3(wave_function_basis(w)),
            size_to_domain(sizes),
        )
        wannier[reset_overflow(k)] = orbitals_from_wave_functions(w, reciprocal_lattice, k)
        w.evc = nothing
    end
    return wannier
end

function wannier_from_unk_dir(
    unk_dir::String,
    wave_functions_list::AbstractVector{WFC},
    domain_scaling_factor::Integer = 2,
)
    k_map, brillouin_zone = i_kpoint_map(wave_functions_list)

    wannier = init_wannier(brillouin_zone)
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

"""
function MMN(mmn_filename::String, ad_hoc_gamma = false)
    file = open(mmn_filename)
    parse_complex(number) =
        parse(Float32, number[1:20]) + parse(Float32, number[21:end]) * 1im
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
    NeighborIntegral(mmn, k_map)

Construct a `NeighborIntegral` from a `MMN` object and a `k_map`.
"""
function NeighborIntegral(mmn::MMN, k_map::Dict{Int64,KPoint})
    n = NeighborIntegral()

    for ((i_kpoint1, i_kpoint2), integral) in mmn.integrals
        kpoint = k_map[i_kpoint1]
        neighbor = add_overflow(k_map[i_kpoint2], mmn.translations[i_kpoint1=>i_kpoint2])
        neighbor = add_overflow(neighbor, -overflow(kpoint))
        kpoint = reset_overflow(kpoint)
        n[kpoint, neighbor] = integral
    end
    return n
end

struct AMN
    n_band::Int64
    n_kpoint::Int64
    n_wannier::Int64
    gauge::Dict{Int64,Matrix{ComplexFxx}}
end

"""
    AMN(amn_filename)

Create an AMN object from a file.
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
            haskey(amn.gauge, i_kpoint) || (amn.gauge[i_kpoint] = zeros(ComplexFxx, (n_band, n_band)))
            amn.gauge[i_kpoint][m, n] = value
        end
    end

    for l_1 = 1:n_kpoint
        Threads.@threads for l_2 = 1:r_1
            line_number = (l_1 - 1) * r_1 + l_2
            m, n, i_kpoint, real_part, complex_part =
                split(strip(lines[line_number]), r"\s+")
            m, n, i_kpoint = (s -> parse(Int64, s)).([m, n, i_kpoint])
            real_part, complex_part = (s -> parse(Float32, s)).([real_part, complex_part])
            update_gauge(i_kpoint, m, n, real_part + complex_part * 1im) 
        end
    end
    close(file)
    return amn
end


"""
    Gauge(grid, amn, k_map, [orthonormalization = true])

Create a gauge from an `AMN` object and a `k_map`.
"""

function Gauge(grid::Grid, amn::AMN, k_map::Dict{Int64,KPoint}, orthonormalization = true)
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
