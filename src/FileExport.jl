export export_to_gaussian_cube

function export_to_gaussian_cube(
    filename::String,
    u::OrbitalSet{UnkBasisOrbital{<:HomeCell}},
    band::AbstractVector{<:Number},
    atoms::Vector = [],
)
    brillouin_zone = grid(u)
    gamma_point = brillouin_zone[0, 0, 0]
    homecell = grid(u[gamma_point][1])

    file = open(filename, "w")
    write(file, "The first line and the second line are comments,\n")
    write(file, "but I don't now what to say.\n")
    origin_x, origin_y, origin_z = cartesian(homecell[0, 0, 0])
    @printf(file, "%i %.8f %.8f %.8f %i\n",
            length(atoms),
            origin_x,
            origin_y,
            origin_z,
            length(band)
        )

    for (s, v) in zip(size(homecell), basis(homecell))
        x, y, z = v.vector
        @printf(file, "%i %.8e %.8e %.8e\n", s, x, y, z)
    end

    function write_atoms(atoms)
        for (atomic_number, charge, x, y, z) in atoms
            @printf(file, "%i %.8e %.8e %.8e %.8e\n", atomic_number, charge, x, y, z)
        end
    end

    isempty(atoms) || write_atoms(atoms)

    function single_wannier_elements(n::Integer)
        raw_elements = zeros(Float32, size(homecell))
        for k in brillouin_zone
            phase = map(r -> exp(1im * (r' * k)), homecell)
            orbital = u[n, k] |> UnkBasisOrbital
            orbital = phase * orbital
            raw_elements += elements(orbital)
        end
        return abs2.(raw_elements)
    end



    wannier_elements = single_wannier_elements.(band)
    n_x, n_y, n_z = size(homecell)

    for m = 1:length(band)
        for i = 1:n_x, j = 1:n_y, k = 1:n_z
            @printf(file, "%.8e ", wannier_elements[m][i, j, k])
            write(file, "\n")
        end
    end

    close(file)
end