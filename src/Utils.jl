export indent, size_to_domain

function size_to_domain(sizes) 
    nx, ny, nz = Int.(sizes)
    return ((-nx ÷ 2 + 1, nx ÷ 2), (-ny ÷ 2 + 1, ny ÷ 2), (-nz ÷ 2 + 1, nz ÷ 2))
end

function indent(str::String)
    lines = split(str, "\n")
    lines = [join(fill(" ", 4)) * line for line in lines]
    return join(lines, "\n")
end
