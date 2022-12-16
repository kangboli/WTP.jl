export Lattice, SimpleCubic, FCC, HCP, make_lattice, enumerate_atoms

abstract type Lattice end

struct SimpleCubic <: Lattice
    a::Float64
    grid::Grid
end

function make_lattice(::Type{SimpleCubic}, a::Number, size::Tuple)
    SimpleCubic(a, make_grid(HomeCell3D, (Vector3(a, 0, 0), Vector3(0, a, 0),
                                          Vector3(0, 0, a)), size))
end

struct FCC <: Lattice
    a::Float64
    grid::Grid
end

function make_lattice(::Type{FCC}, a::Number, size::Tuple)
    FCC(a, make_grid(HomeCell3D, (Vector3(0, a/2, a/2), Vector3(a/2, 0, a/2),
                                  Vector3(a/2, a/2, 0)), size))
end

struct HCP <: Lattice
    a::Float64
    c::Float64
    grid::Grid
end

function make_lattice(::Type{HCP}, a::Number, c::Number, size::Tuple)
    HCP(a, c, make_grid(HomeCell3D, (Vector3(a, 0, 0), Vector3(1/2 * a, sqrt(3)/2 * a, 0), Vector3(0, 0, c)), size))
end


grid(l::Lattice) = l.grid


function enumerate_atoms(lattice::HCP, atomic_number::Int)
    origin_points = coefficients.(collect(grid(lattice)))
    
    function create_atom_pair(p)
    a = make_atom(atomic_number, grid(lattice), p) 
    b = make_atom(atomic_number, grid(lattice), [p[1]+1/3, p[2]+1/3,
                                                           p[3]+1/2])
        return [a, b]
    end

    return vcat(create_atom_pair.(origin_points)...)
end


function enumerate_atoms(lattice::Lattice, atomic_number::Int)
    origin_points = coefficients.(collect(grid(lattice)))
    return [make_atom(atomic_number, grid(lattice), p) for p in origin_points]
end

