using FortranFiles

export UNK,
    WFC,
    MMN,
    AMN,
    load_evc!,
    single_orbital_from_wave_functions,
    orbitals_from_wave_functions,
    brillouin_zone_from_k_coords,
    WannierFromSave,
    wave_functions_from_directory,
    estimate_sizes,
    wave_functions_basis,
    wannier_from_save,
    orbitals_from_unk,
    wannier_from_unk_dir,
    i_kpoint_map

"""
Ugly code for dealing with ugly files.
"""
const parse_line = (line, type) -> (n -> parse(type, n)).(split(strip(line), r"\s+"))

# This factor is possibly either from Nyquist-Shannon sampling theorem.
# or from the energy cutoff for the density. 
# Settting this to 2 gives the same FFT as QE.
const DOMAIN_SCALING_FACTOR = 2

struct UNK
    nx::Int64
    ny::Int64
    nz::Int64
    k::Int64
    n_band::Int64
    psi_r::Matrix{ComplexF64}
    filename::String
end



"""
The wave functions are read such that each row is an orbital
to avoid transposing in memory for SCDM.


This takes 10s to load a 1.4G unk file using 12 threads and an NVMe drive.
"""
function UNK(unk_filename::String)
    # Read the header information.
    file = open(unk_filename)
    parse_complex(number) =
        parse(Float32, number[1:20]) + parse(Float32, number[21:40]) * 1im
    nx, ny, nz, k, n_band = parse_line(readline(file), Int64)
    unk = UNK(nx, ny, nz, k, n_band, zeros(ComplexF64, n_band, nx * ny * nz), unk_filename)
    ng = nx * ny * nz
    lines = readlines(file)
    for n = 1:n_band
        Threads.@threads for i = 1:ng
            unk.psi_r[n, i] = parse_complex(lines[ng*(n-1)+i])
        end
    end
    close(file)
    return unk
end

mutable struct WFC
    filename::String
    i_kpoint::Int64
    k_coords::Vector{Float64}
    i_spin::Int64
    gamma_only::Bool
    scale_factor::Float64
    n_planewaves::Int64
    max_n_planewaves::Int64
    n_polerizations::Int64
    n_band::Int64

    b1::Vector{Float64}
    b2::Vector{Float64}
    b3::Vector{Float64}
    miller::Matrix{Int64}
    evc::Union{Matrix{ComplexF64},Nothing}
end

wave_functions_basis(wave_functions::WFC) = [wave_functions.b1 wave_functions.b2 wave_functions.b3]

"""
This may not work when Quantum Espresso is compiled with a different Fortran
compiler.  Fortran does not enforce a standard on the size of it's primitive
types in binary files.  The ways to get consistent behaviors are

1. Stay away from Fortran binary files. 
"""
function WFC(wave_function_filename::String)
    # See the format in: https://gitlab.com/QEF/q-e/-/wikis/Developers/Format-of-wave_functions-files 
    data = FortranFile(wave_function_filename)
    i_kpoint, k_coords, i_spin, gamma_only, scale_factor =
        read(data, Int32, (Float64, 3), Int32, Int32, Float64)
    gamma_only = gamma_only > 0
    n_planewaves, max_n_planewaves, n_polerizations, n_band =
        read(data, Int32, Int32, Int32, Int32)
    b1, b2, b3 = read(data, (Float64, 3), (Float64, 3), (Float64, 3))
    miller = read(data, (Int32, (3, max_n_planewaves)))
    close(data)

    return WFC(
        wave_function_filename,
        i_kpoint,
        k_coords,
        i_spin,
        gamma_only,
        scale_factor,
        n_planewaves,
        max_n_planewaves,
        n_polerizations,
        n_band,
        b1,
        b2,
        b3,
        miller,
        nothing,
    )

end

"""
Load the evc (frequency space wave-functions) into wave_functions.
"""
function load_evc!(w::WFC)
    data = FortranFile(w.filename)
    for _ = 1:4
        read(data)
    end
    w.evc = zeros(ComplexF64, w.n_polerizations * w.max_n_planewaves, w.n_band)
    for i = 1:w.n_band
        w.evc[:, i] = read(data, (ComplexF64, w.n_polerizations * w.max_n_planewaves))
    end
    close(data)
end

function i_kpoint_map(wave_functions_list::AbstractVector{WFC})
    kmap = Dict{Int64, KPoint}()

    c = (w->w.k_coords).(wave_functions_list)
    brillouin_zone = brillouin_zone_from_k_coords(c,
        wave_functions_basis(wave_functions_list[1]))
    for w in wave_functions_list
        kmap[w.i_kpoint] = snap(brillouin_zone, w.k_coords)
    end
    return kmap, brillouin_zone
end

"""
Create an unk orbital from a WFC structure.

k can be out of the first Brillouin zone. An image withing the first 
Brillouin zone will be constructed and returned.
"""
function single_orbital_from_wave_functions(w::WFC, l::ReciprocalLattice, k::KPoint, n::Int)
    w.evc !== nothing || error("The wave functions are not loaded.")
    elements = zeros(ComplexF64, size(l)...)
    orbital = UnkBasisOrbital(l, elements, reset_overflow(k), n)

    i_kpoint!(orbital, w.i_kpoint)

    for i = 1:w.max_n_planewaves
        g = grid_vector_constructor(l, w.miller[:, i] + overflow(k))
        orbital[g] = w.evc[i, n]
    end

    return orbital
end

"""
Create an unk orbital from a UNK structure.

k can be out of the first Brillouin zone. The resulting orbital 
will also be an image outside the first Brillouin zone.
"""
function single_orbital_from_unk(unk::UNK, homecell::HomeCell, k::KPoint, index_band::Int)
    # TODO: Figure out the layout.
    if (unk.nx != unk.ny) || (unk.nx != unk.nz)
        @warn "UNK files with different nx ny nz dimension has not been
                implemented because pw2wannier90.x has no documentation"
    end

    elements = reshape(unk.psi_r[index_band, :], (unk.nx, unk.ny, unk.nz))
    UnkBasisOrbital(homecell, elements, k, index_band)
end

"""
Create all the orbitals that are in the wave_functions structure.

k can be out of the first Brillouin zone. Images within the first Brillouin zone
will be constructed and returned.
"""
function orbitals_from_wave_functions(wave_functions::WFC, reciprocal_lattice::ReciprocalLattice, k::KPoint)
    return [
        single_orbital_from_wave_functions(wave_functions, reciprocal_lattice, k, index_band) for
        index_band = 1:wave_functions.n_band
    ]
end

"""
The returned orbital will be in the first Brillouin zone.
"""
function orbitals_from_unk(unk::UNK, homecell::HomeCell, k::KPoint)

    """
    Convert the orbital from its image in Wannier brillouin zone to QE brillouin zone.
    """
    function standard_brillouin_zone(orbital, k::KPoint)
        has_overflow(k) || return orbital
        return standardize(translate(orbital, grid(orbital)[overflow(k)...]))
    end

    function reciprocal_orbital(b)
        o = fft(single_orbital_from_unk(unk, homecell, k, b))
        wtp_normalize!(o)
        o = standard_brillouin_zone(o, k)
    end
    reciprocal_orbital.(collect(1:unk.n_band))
end

estimate_sizes(wave_functions::WFC, i::Int) =
    2 * DOMAIN_SCALING_FACTOR * (max(abs.(wave_functions.miller[i, :])...))

"""
Construct a BrillouinZone (grid) from a set of kpoints and the reciprocal basis vectors.

This is done by first writing the kpoints in terms of the reciprocal basis.
Then we take the three smallest kpoints besides the gamma point. 
"""
function brillouin_zone_from_k_coords(
    k_coords::Vector{Vector{Float64}},
    reciprocal_basis::Matrix{Float64},
)
    gamma_point = k_coords[argmin(norm.(k_coords))]
    offset_kpoints = [k_coords[i] - gamma_point for i = 1:length(k_coords)]
    # Represents the kpoint in the reciprocal basis.
    in_reciprocal_basis = [reciprocal_basis \ k for k in offset_kpoints]

    find_positive_min(cs::Vector{Float64}) =
        min(filter((c) -> !isapprox(c, 0, atol = 1e-8), cs)...)

    # The points right next to the gamma point gives the size of the grid.
    brillouin_sizes = Tuple(
        round(1 / find_positive_min([kr[i] for kr in in_reciprocal_basis])) for i = 1:3
    )
    BrillouinZone(
        matrix_to_vector3(reciprocal_basis * inv(diagm([brillouin_sizes...]))),
        size_to_domain(brillouin_sizes),
    )
end

function wave_functions_from_directory(save_dir::String)
    is_wave_functions(f::String) = startswith(basename(f), "wfc") && endswith(basename(f), ".dat")
    return [WFC(file) for file in readdir(save_dir, join = true) if is_wave_functions(file)]
end


function wannier_from_save(wave_functions_list::AbstractVector{WFC})
    # wave_functions_list = wave_functions_from_directory(save_dir)
    k_map, brillouin_zone = i_kpoint_map(wave_functions_list)

    # c = (w->w.k_coords).(wave_functions_list)
    # brillouin_zone =
    #     brillouin_zone_from_k_coords(c, wave_functions_basis(wave_functions_list[1]))
    wannier = Wannier(brillouin_zone)
    sizes = Tuple(maximum((w) -> estimate_sizes(w, i), wave_functions_list) for i = 1:3)

    for w in wave_functions_list
        k = k_map[w.i_kpoint]
        # elements(gauge(wannier))[miller_to_standard(k, translation(wannier))...] =
        gauge(wannier)[reset_overflow(k)] = Matrix{Float64}(I, w.n_band, w.n_band)

        load_evc!(w)
        reciprocal_lattice = ReciprocalLattice(matrix_to_vector3(wave_functions_basis(w)), size_to_domain(sizes))
        # elements(wannier)[miller_to_standard(k, translation(wannier))...] =
        wannier[reset_overflow(k)] = orbitals_from_wave_functions(w, reciprocal_lattice, k,)
    end
    return wannier
end

function wannier_from_unk_dir(unk_dir::String, wave_functions_list::AbstractVector{WFC})
    k_map, brillouin_zone = i_kpoint_map(wave_functions_list)

    wannier = Wannier(brillouin_zone)
    sizes = Tuple(maximum((w) -> estimate_sizes(w, i), wave_functions_list) for i = 1:3)

    for w in wave_functions_list
        k = k_map[w.i_kpoint]
        unk = UNK("$(unk_dir)/UNK$(lpad(w.i_kpoint, 5, "0")).1")
        # wannier.gauge[miller_to_standard(k, translation(wannier))...] =
        gauge(wannier)[reset_overflow(k)] = Matrix{Float64}(I, unk.n_band, unk.n_band)

        # println("reading k point: $(coefficients(k))\n ik: $(w.i_kpoint)")
        reciprocal_lattice =
            ReciprocalLattice(matrix_to_vector3(wave_functions_basis(w)), size_to_domain(sizes))
        homecell = transform_grid(reciprocal_lattice)
        reciprocal_orbitals = orbitals_from_unk(unk, homecell, k)
        (o -> i_kpoint!(o, w.i_kpoint)).(reciprocal_orbitals)
        # wannier.elements[miller_to_standard(k, translation(wannier))...] =
        wannier[reset_overflow(k)] = reciprocal_orbitals

    end

    return wannier
end

struct MMN
    n_band::Int64
    n_kpoint::Int64
    n_neighbor::Int64
    integrals::Dict{Pair{Int, Int}, Matrix{ComplexF64}}
    translations::Dict{Pair{Int, Int}, Tuple{Int, Int, Int}}
end

"""
From Wannier90 user_guide.pdf:

    The last three integers specify the G vector, in reciprocal lattice units,
    that brings the k-point specified by the second integer, and that thus lives inside
    the first BZ zone, to the actual k + b we need.
"""
function MMN(mmn_filename::String)
    file = open(mmn_filename)
    parse_complex(number) =
        parse(Float32, number[1:20]) + parse(Float32, number[21:end]) * 1im
    _ = readline(file) # Skip the date.
    # println(date)
    n_band, n_kpoint, n_neighbor = parse_line(readline(file), Int64)
    mmn = MMN(n_band, n_kpoint, n_neighbor, Dict(), Dict())
    lines = readlines(file)

    r1, r2 = n_neighbor*(n_band^2+1), (n_band^2+1)
    for l1 = 1:n_kpoint
        # Threads.@threads 
        for l2 = 1:n_neighbor
            start = (l1-1) * r1 + (l2-1) * r2 + 1
            i_kpoint, i_neighbor, gx, gy, gz = parse_line(lines[start], Int64)
            mmn.translations[i_kpoint=>i_neighbor] = (gx, gy, gz)
            mmn.integrals[i_kpoint=>i_neighbor] = 
            [parse_complex(lines[start + (i-1)*n_band+j]) for j=1:n_band, i=1:n_band]
        end
    end
    close(file)
    return mmn
end

"""
"""
function NeighborIntegral(mmn::MMN, k_map::Dict{Int64, KPoint})
    n = NeighborIntegral()

    for ((i_kpoint1, i_kpoint2), integral) in mmn.integrals
        kpoint = k_map[i_kpoint1]
        neighbor = add_overflow(k_map[i_kpoint2],  mmn.translations[i_kpoint1=>i_kpoint2])
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
    gauge::Dict{Int64, Matrix{ComplexF64}}
end

function AMN(amn_filename::String)
    file = open(amn_filename)
    _ = readline(file) # Skip the date.
    n_band, n_kpoint, n_wannier = parse_line(readline(file), Int64)
    amn = AMN(n_band, n_kpoint, n_wannier, Dict())
    lines = readlines(file)
    r_1 = n_band^2
    # Threads.@threads 
    for l_1 = 1:n_kpoint
        for l_2 = 1:r_1
            line_number = (l_1-1)*r_1+l_2
            m, n, i_kpoint, real_part, complex_part = split(strip(lines[line_number]), r"\s+")
            m, n, i_kpoint = (s->parse(Int64, s)).([m, n, i_kpoint])
            real_part, complex_part = (s->parse(Float32, s)).([real_part, complex_part])
            haskey(amn.gauge, i_kpoint) || (amn.gauge[i_kpoint] = zeros(ComplexF64, (n_band, n_band)))
            amn.gauge[i_kpoint][m, n] = real_part + complex_part * 1im
        end
    end
    close(file)
    return amn
end

function gauge(amn::AMN, k_map::Dict{Int64, KPoint})
    g = Dict{KPoint, Matrix{ComplexF64}}()

    for (i_kpoint, gauge) in amn.gauge
        g[k_map[i_kpoint]] = gauge
    end
    return g
end
    # file = open(unk_filename)
    # nx, ny, nz, k, n_band = parse_line(readline(file), Int64)
    # unk = UNK(nx, ny, nz, k, n_band, zeros(ComplexF64, n_band, nx * ny * nz), unk_filename)
    # ng = nx * ny * nz
    # lines = readlines(unk_filename)
    # for n = 1:n_band
    #     Threads.@threads for i = 1:ng
    #         unk.psi_r[n, i] = parse_complex(lines[ng*(n-1)+i+1])
    #     end
    # end
    # close(file)