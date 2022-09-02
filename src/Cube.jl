export Cube, to_cube, add_atom!, write_cube

"""
Data structure for dealing with the cube file
"""
abstract type CubeUnit end
struct CubeAngstrom <: CubeUnit end
struct CubeBohr <: CubeUnit end

mutable struct Cube
    # The first two lines of the header are comments, they are generally ignored by
    # parsing packages or used as two default labels.
    header::Vector{String}
    # The third line has the number of atoms included in the file followed by the
    # position of the origin of the volumetric data.
    n_atom::Int
    origin::Vector{Float64}
    # The next three lines give the number of voxels along each axis (x, y, z)
    # followed by the axis vector. Note this means the volume need not be aligned
    # with the coordinate axis, indeed it also means it may be sheared although
    # most volumetric packages won't support that. The length of each vector is the
    # length of the side of the voxel thus allowing non cubic volumes. If the sign
    # of the number of voxels in a dimension is positive then the units are Bohr,
    # if negative then Angstroms.
    axis::Vector{Tuple{CubeUnit, UInt, Vector{Float64}}}
    # The last section in the header is one line for each atom consisting of 5
    # numbers, the first is the atom number, the second is the charge, and the
    # last three are the x,y,z coordinates of the atom center.
    atoms::Vector{Tuple{UInt, Float64, Vector{Float64}}}
    # Volumetric data The volumetric data is straightforward, one floating
    # point number for each volumetric element. The original Gaussian format
    # arranged the values in the format shown below in the example, most
    # parsing programs can read any white space separated format. Traditionally
    # the grid is arranged with the x axis as the outer loop and the z axis as
    # the inner loop, for example, written as
    volumn::Array{Float64, 3}
end


function to_cube(on_grid::OnGrid, atoms=[])
    g = grid(on_grid)
    c = Cube(
        ["The first two lines must tell a joke.\n", 
         "The parser won't work unless tickled.\n"],
        0,
        [0.0, 0.0, 0.0],
        map((s, b) -> (CubeBohr(), s, b), size(g), eachcol(basis_matrix(g))),
        Vector{Tuple{UInt,Float64,Vector{Float64}}}(),
        elements(on_grid),
    )
    map(a->add_atom!(c, a),  atoms)
    
    return c
end


function add_atom!(c::Cube, atom::AbstractAtom)
    c.n_atom += 1
    push!(c.atoms, (atomic_number(atom), get(atom.meta, "charge", 0.0), position(atom)))
end

function write_cube(c::Cube, filename::String)
    f = open(filename, "w")
    map(line -> write(f, line), c.header)
    @printf(f, "%i %.8f %.8f %.8f\n", c.n_atom, c.origin...)
    unit_sign(u::CubeUnit) = u == CubeBohr() ? "" : "-"
    map(ax -> @printf(f, "%s%i %.8e %.8e %.8e\n", unit_sign(ax[1]), ax[2], ax[3]...), c.axis)
    map(a -> @printf(f, "%i %.8e %.8e %.8e %.8e\n", a[1], a[2], a[3]...), c.atoms)
    join(map(v -> @sprintf("%.8e\n", v), c.volumn)) |> s -> write(f, s)
    write(f, "\n")
    close(f)
end
