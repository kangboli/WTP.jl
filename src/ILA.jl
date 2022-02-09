
export ApproximationScheme,
    CosScheme,
    CosScheme3D,
    finite_difference,
    find_shells,
    compute_weights,
    spread,
    shells,
    weights,
    populate_integral_table!,
    center,
    second_moment,
    spread,
    neighbor_basis_integral,
    BranchStable,
    W90BranchCut,
    TruncatedConvolution,
    gauge_gradient,
    ILAOptimizer,
    logger,
    scheme,
    NeighborIntegral,
    integrals,
    FletcherReeves

"""
The neighbor integrals indexed by two kpoints.  For each pair of kpoint, the
integrals are stored as a matrix. This is the same as the ``M_{mn}^{k, b}``
matrix in MLWF. The matrix elements are accessed by `M[k, k+b][m, n]`
"""
struct NeighborIntegral
    integrals::Dict{Pair{KPoint,KPoint},Matrix{ComplexFxx}}
end

abstract type ApproximationScheme end

abstract type CosScheme <: ApproximationScheme end

struct CosScheme3D <: CosScheme
    neighbor_shells::AbstractVector{Vector{<:KPoint}}
    weights::AbstractVector{Number}
    neighbor_basis_integral::NeighborIntegral
end

"""
    CosScheme3D(ũ, n_shells = 1)

The function ``r^2`` can be approximated with ``w_{\\mathbf{b}} \\cos(\\mathbf{b}^T
\\mathbf{r})`` functions, which are then used for approximating the convolution
between ``r^2`` and `u` (inverse fft of `ũ`).
The `CosScheme3D` includes shells of `\\mathbf{b}` vectors and their 
corresponding weights ``w_{\\mathbf{b}}``.

Example: 

```jldoctest orbital_set
julia> scheme = CosScheme3D(ũ);

julia> length(shells(scheme))
1
```
"""
function CosScheme3D(u::OrbitalSet{UnkBasisOrbital{ReciprocalLattice3D}}, n_shells = 1)
    neighbor_shells = find_shells(grid(u), n_shells)
    weights = compute_weights(neighbor_shells)
    weights === nothing && return CosScheme3D(u, n_shells + 1)
    neighbor_integral = NeighborIntegral()

    scheme = CosScheme3D(neighbor_shells, weights, neighbor_integral)
    populate_integral_table!(scheme, u)
    return scheme
end


"""
    shells(scheme)

Shells of ``\\mathbf{b}`` vectors involved in the approximation scheme.  Each
shell is a vector of kpoints.

Example: 

```jldoctest orbital_set
julia> shells(scheme)

1-element Vector{Vector{<:KPoint}}:
 KPoint[GridVector{BrillouinZone3D}:
    coefficients: [-1, -1, -1]
, GridVector{BrillouinZone3D}:
    coefficients: [-1, 0, 0]
, GridVector{BrillouinZone3D}:
    coefficients: [0, -1, 0]
, GridVector{BrillouinZone3D}:
    coefficients: [0, 0, -1]
, GridVector{BrillouinZone3D}:
    coefficients: [0, 0, 1]
, GridVector{BrillouinZone3D}:
    coefficients: [0, 1, 0]
, GridVector{BrillouinZone3D}:
    coefficients: [1, 0, 0]
, GridVector{BrillouinZone3D}:
    coefficients: [1, 1, 1]
]
```
"""
shells(scheme::ApproximationScheme)::AbstractVector{Vector{<:KPoint}} =
    scheme.neighbor_shells

"""
    weights(scheme)

The weights corresponding to each shell within a scheme. The weights are ordered
from the inner-most to the outer-most shell.

```jldoctest orbital_set
julia> weights(scheme)
1-element Vector{Number}:
 5.336038037571918
```
"""
weights(scheme::ApproximationScheme) = scheme.weights

"""
    neighbor_basis_integral(scheme)

The integrals between neighboring k-points (The MMN matrix). 
The integral is amonst immediate neighbor because the ``\\cos`` approximation
is truncated at the first mode.

```jldoctest orbital_set
julia> M = neighbor_basis_integral(scheme)
julia> M[brillouin_zone[0, 0, 0], brillouin_zone[0, 0, 1]]
4×4 Matrix{ComplexF64}:
     0.85246+0.512198im     0.0798774-0.00342381im  4.54441e-8-1.98313e-8im  8.32139e-9+4.06675e-8im
   0.0057819+0.021992im     -0.363241-0.324958im     -0.393272-0.187671im     -0.489628+0.127022im
   0.0112149-0.000481798im  -0.194858+0.141126im      0.665299+0.0837434im    -0.498015-0.0113525im
 -0.00992131-0.0216156im     0.432817+0.269315im     -0.319281+0.193723im     -0.524017-0.0193465im
```
"""
neighbor_basis_integral(scheme::ApproximationScheme) = scheme.neighbor_basis_integral
# transformed_integral(scheme::ApproximationScheme) = scheme.transformed_integral

"""
    find_neighbors(k, scheme)

Find the set of relevant neighbors for a kpoint under a scheme.
"""
function find_neighbors(kpoint::KPoint, scheme::ApproximationScheme)
    dk_list = vcat(shells(scheme)...)
    # TODO: The negative part may not be necessary.
    return (dk -> kpoint + dk).([dk_list; -dk_list])
end

# function gauge_transform(scheme::ApproximationScheme, U::Gauge)
#     @set scheme.neighbor_basis_integral = gauge_transform(neighbor_basis_integral(scheme), U)
# end


"""
    NeighborIntegral() 

The neighbor integrals is roughly a dictionary with pairs of neighboring
k-points as the keys. One can create one just by
"""
NeighborIntegral() = NeighborIntegral(Dict())
integrals(n::NeighborIntegral) = n.integrals

function Base.hash(p::Pair{<:KPoint,<:KPoint})
    m, n = p
    a_nice_prime_number = 7
    return hash(m) * a_nice_prime_number + hash(n)
end

function Base.:(==)(p_1::Pair{<:KPoint,<:KPoint}, p_2::Pair{<:KPoint,<:KPoint})
    m_1, n_1 = p_1
    m_2, n_2 = p_2
    return m_1 == m_2 && n_1 == n_2
end

"""
    getindex(M, k_1, k_2)

The integral matrix between the neighboring k-points `k_1` and `k_2`.
Can also write `M[k_1, k_2]`. Note that `M[k_1, k_2] = M[k_2, k_1]'`.
So only one of the matrices is stored.

```jldoctest orbitla_set
julia> M[brillouin_zone[0, 0, 0], brillouin_zone[0, 0, -1]] == M[brillouin_zone[0, 0, -1], brillouin_zone[0, 0, 0]]'
true
```
"""
function Base.getindex(neighbor_integral::NeighborIntegral, k_1::KPoint, k_2::KPoint)
    # coefficients(k_1) == coefficients(k_2) && return I
    i = integrals(neighbor_integral)
    # haskey(i, k_1 => k_2) && return i[k_1=>k_2]
    # return adjoint(i[k_2=>k_1])
    result = get(i, k_1 => k_2, nothing)
    return result === nothing ? adjoint(i[k_2=>k_1]) : result
end

"""
    setindex!(M, value, g...)

Set the integral matrix between two k-points `g[1]` and `g[2]` to
`value`. Can also write `M[g...] = value`.
"""
function Base.setindex!(neighbor_integral::NeighborIntegral, value::AbstractMatrix, g::Vararg{<:KPoint})
    g_1, g_2 = g
    integrals(neighbor_integral)[g_1=>g_2] = value
end

"""
    haskey(M, k_1, k_2)

Check if the integral matrix between `k_1` and `k_2` has been computed and stored.

```jldoctest orbital_set
julia> haskey(M, brillouin_zone[0, 0, 1], brillouin_zone[0, 0, 0])
true

julia> haskey(M, brillouin_zone[0, 0, 1], brillouin_zone[0, 0, -1])
false
```
"""
function Base.haskey(neighbor_integral::NeighborIntegral, k_1::KPoint, k_2::KPoint)
    i = integrals(neighbor_integral)
    return haskey(i, k_1 => k_2) || haskey(i, k_2 => k_1)
end

"""
    gauge_transform(M, gauge)

Perform a gauge transform on the neighbor integrals.

``U^{k \\dagger} M^{k, k+b} U^{k+b}``

```jldoctest orbital_set
julia> M = gauge_transform(M, U);
julia> M[brillouin_zone[0, 0, 1], brillouin_zone[0, 0, 0]]
4×4 adjoint(::Matrix{ComplexF64}) with eltype ComplexF64:
    0.85246-0.512198im    0.0057819-0.021992im  0.0112149+0.000481798im  -0.00992131+0.0216156im
  0.0798774+0.00342381im  -0.363241+0.324958im  -0.194858-0.141126im        0.432817-0.269315im
 4.54441e-8+1.98313e-8im  -0.393272+0.187671im   0.665299-0.0837434im      -0.319281-0.193723im
 8.32139e-9-4.06675e-8im  -0.489628-0.127022im  -0.498015+0.0113525im      -0.524017+0.0193465im
```

"""
function gauge_transform(neighbor_integral::NeighborIntegral, gauge::Gauge)
    t = NeighborIntegral()
    for ((k_1, k_2), integral) in integrals(neighbor_integral)
        t[k_1, k_2] = adjoint(gauge[k_1]) * integral * gauge[k_2]
    end
    return t
end


"""
The Brillouin zone on which the finite difference scheme is defined.
"""
grid(scheme::ApproximationScheme) = grid(shells(scheme)[1][1])

function find_shells(grid::Grid, n_shell::Int)
    shells = SortedDict{Real,Vector{KPoint}}()
    d = collect(-2*n_shell:2*n_shell)
    for i in d, j in d, k in d
        k = make_grid_vector(grid, [i, j, k])
        key = round(norm(cartesian(k)), digits = 5)
        haskey(shells, key) ? append!(shells[key], [k]) : shells[key] = [k]
    end

    return collect(values(shells))[2:n_shell+1]
end

"""
Solve Aw = q
"""
function compute_weights(neighbor_shells::Vector{Vector{T}}) where {T<:AbstractGridVector}

    indices = SortedDict(
        (1, 1) => 1,
        (1, 2) => 2,
        (1, 3) => 3,
        (2, 2) => 4,
        (2, 3) => 5,
        (3, 3) => 6,
    )

    A = zeros((6, length(neighbor_shells)))
    c(b, i) = cartesian(b)[i]
    for s = 1:length(neighbor_shells)
        A[:, s] =
            [sum((b -> c(b, i) * c(b, j)).(neighbor_shells[s])) for (i, j) in keys(indices)]
    end

    q = [i == b ? 1 : 0 for (i, b) in keys(indices)]
    w = A \ q

    return isapprox(A * w, q, atol = 1e-5) ? w : nothing
end

"""
Use the same algorithm as in the original Wannier90 (without guiding centers,
fixed centers, etc.).  Apply this algorithm if you like to cross validate with
Wannier90.
"""
abstract type W90BranchCut end
"""
Mostly the same algorithms as `W90BranchCut`, but make a potential more
consistent choice of branch cut.
"""
abstract type BranchStable end
"""
The truncated convolution algorithm. We believe this to be the state of the art as of 2022.
"""
abstract type TruncatedConvolution end

center(scheme::CosScheme, n::Integer) = center(scheme, n, W90BranchCut)
center(scheme::CosScheme, n::Integer, ::Type{T}) where {T} = center(neighbor_basis_integral(scheme), scheme, n, T)
center(M::NeighborIntegral, scheme::CosScheme, n::Integer) = center(M, scheme, n, W90BranchCut)

"""
    center(M, scheme, n, ::Type{W90BranchCut})

Compute the center of the `n`th Wannier orbital using the original Wannier90 approach.
The replication is exact.

Example:

```jldoctest orbital_set
julia> center(M, scheme, 1, W90BranchCut)
3-element Vector{Float64}:
 -8.73495454038011
  3.9151488754936588
  4.021344063664258
julia> center(M, scheme, 2, W90BranchCut)
3-element Vector{Float64}:
  0.11705358388655285
 -1.4773864607361118
 -2.4684488265569624
```
"""
function center(M::NeighborIntegral, scheme::CosScheme, n::Int, ::Type{W90BranchCut})
    kpoint_contribution(k::KPoint) = -sum(zip(weights(scheme), shells(scheme))) do (w, shell)
        sum(b -> w * cartesian(b) * angle(M[k, k+b][n, n]), shell)
    end

    brillouin_zone = collect(grid(scheme))
    return sum(kpoint_contribution.(brillouin_zone)) / length(brillouin_zone)
end

"""
    center(M, scheme, n, ::Type{TruncatedConvolution})

Compute the center of the `n`th Wannier orbital using the turncated convolution algorihtm.

Example:

```jldoctest orbital_set
julia> center(M, scheme, 1, TruncatedConvolution)
3-element Vector{Float64}:
 -8.70234629129622
  4.016099860574682
  4.093053970246974
julia> center(M, scheme, 1, TruncatedConvolution)
3-element Vector{Float64}:
  1.3634042350006577
 -2.712474285697622
 -3.6252083684030336
```

One can compare this with the "exact" center

```jldoctest orbital_set
julia> wanniers = commit_gauge(ũ)(:);
julia> _, r̃2 = compute_r2(supercell(u));
julia> c_1, σ²_1 = center_spread(fft(abs2(ifft(wanniers[1])), false), r̃2)
([-8.831796622851812, 3.953063696051289, 3.940388714312329], 28.748708554079474)
julia> c_2, σ²_2 = center_spread(fft(abs2(ifft(wanniers[2])), false), r̃2)
([-9.192231000846695, 7.604578700842824, 6.189389913055544], 27.460489592320883)
```
"""
function center(M::NeighborIntegral, scheme::CosScheme, n::Int, ::Type{TruncatedConvolution})
    brillouin_zone = collect(grid(scheme))

    phase(b::KPoint) = angle(1 / length(brillouin_zone) * sum(k -> M[k, k+b][n, n], brillouin_zone))

    return -sum(zip(weights(scheme), shells(scheme))) do (w, shell)
        sum(unique(k -> Set([k, -k]), shell)) do b
            ϕ⁺, ϕ⁻ = phase(b), phase(-b)
            branch = (sign(ϕ⁺) == sign(ϕ⁻) ? -1 : 1)
            w * cartesian(b) * ϕ⁺ + w * cartesian(-b) * branch * ϕ⁻
        end
    end
end


function center(M::NeighborIntegral, scheme::CosScheme, n::Int, ::Type{BranchStable})
    kpoint_contribution(k::KPoint) = -sum(zip(weights(scheme), shells(scheme))) do (w, shell)
        sum(unique(k -> Set([k, -k]), shell)) do b
            ϕ⁺ = M[k, k+b][n, n] |> angle
            ϕ⁻ = M[k, k-b][n, n] |> angle
            branch = (sign(ϕ⁺) == sign(ϕ⁻) ? -1 : 1)
            w * cartesian(b) * ϕ⁺ + w * cartesian(-b) * branch * ϕ⁻
        end
    end

    brillouin_zone = scheme |> grid |> collect
    return sum(kpoint_contribution.(brillouin_zone)) / prod(size(brillouin_zone))

end

second_moment(scheme::CosScheme, n::Int) = second_moment(neighbor_basis_integral(scheme), scheme, n)

function second_moment(M::NeighborIntegral, scheme::CosScheme, n::Int)
    kpoint_contribution(k::KPoint) =
        sum(zip(weights(scheme), shells(scheme))) do (w, shell)
            sum(b -> w * (1 - abs2(M[k, k+b][n, n]) + angle(M[k, k+b][n, n])^2), shell)
        end
    brillouin_zone = collect(grid(scheme))
    return sum(kpoint_contribution.(brillouin_zone)) / prod(size(brillouin_zone))
end

spread(scheme::CosScheme, n::Integer) = spread(scheme, n, W90BranchCut)
spread(scheme::CosScheme, n::Integer, ::Type{T}) where {T} = spread(neighbor_basis_integral(scheme), scheme, n, T)
spread(M::NeighborIntegral, scheme::CosScheme, n::Integer) = spread(M, scheme, n, W90BranchCut)

"""
    spread(M, scheme, n, W90BranchCut)

Compute the spread of the `n`th Wannier orbital using the original Wannier90 approach.

```jldoctest orbitla_set
julia> spread(M, scheme, 1, W90BranchCut)
15.742034477681969
julia> spread(M, scheme, 2, W90BranchCut) # failure.
133.09413338071354
```
"""
spread(M::NeighborIntegral, scheme::CosScheme, n::Integer, ::Type{T}) where {T} =
    second_moment(M, scheme, n) - norm(center(M, scheme, n, T))^2

"""
    spread(M, scheme, n, TruncatedConvolution)

Compute the spread of the `n`th Wannier orbital using the truncated convolution.

```jldoctest orbitla_set
julia> spread(M, scheme, 1, TruncatedConvolution)
17.50438313709964
julia> spread(M, scheme, 2, TruncatedConvolution)
17.313972338201154
```

Compare this to the "exact" spread.

```jldoctest orbital_set 
julia> σ²_1, σ²_2
(28.748708554079474, 27.460489592320883)
```
"""
function spread(M::NeighborIntegral, scheme::CosScheme, n::Integer, ::Type{TruncatedConvolution})
    brillouin_zone = collect(grid(scheme))
    ρ̃(b) = sum(k -> M[k, k+b][n, n], brillouin_zone) / length(brillouin_zone)

    sum(zip(weights(scheme), shells(scheme))) do (w, shell)
        sum(b -> 2w * (1 - abs(ρ̃(b))), shell)
    end
end


function all_spread(U, scheme, ::Type{TruncatedConvolution})
    M = gauge_transform(neighbor_basis_integral(scheme), U)
    N = n_band(M)
    brillouin_zone = collect(grid(scheme))
    ρ̃(b) = sum(k -> diag(M[k, k+b]), brillouin_zone) / length(brillouin_zone)

    sum(zip(weights(scheme), shells(scheme))) do (w, shell)
        sum(b -> 2w * (ones(N) - abs.(ρ̃(b))), shell)
    end
end

function populate_integral_table!(scheme::ApproximationScheme, u::OrbitalSet)
    brillouin_zone = grid(u)
    M = neighbor_basis_integral(scheme)

    """
    Compute the mmn matrix between k1 and k2.
    """
    function mmn_matrix(k_1::T, k_2::T) where {T<:AbstractGridVector{<:BrillouinZone}}
        U = hcat(vectorize.(u[k_1])...)
        V = hcat(vectorize.(u[k_2])...)
        return adjoint(U) * V
        # return [braket(m, n) for m in u[k_1], n in u[k_2]]
    end

    @showprogress for k in collect(brillouin_zone)
        for neighbor in find_neighbors(k, scheme)
            # M[k, neighbor] = adjoint(U[k]) * mmn_matrix(k, neighbor) * U[neighbor]
            haskey(M, k, neighbor) && continue
            M[k, neighbor] = mmn_matrix(k, neighbor)
        end
    end
    return M
end

# for k in brillouin_zone
#     (m -> cache!(m, M)).(u[k])
# end

function n_band(M::NeighborIntegral)
    first_matrix = collect(values(integrals(M)))[1]
    return size(first_matrix, 1)
end

gauge_gradient(U::Gauge, scheme::CosScheme, brillouin_zone::BrillouinZone) =
    gauge_gradient(U::Gauge, scheme::CosScheme, brillouin_zone::BrillouinZone, W90BranchCut)


"""
    gauge_gradient(U, scheme, brillouin_zone, W90BranchCut)

Gauge gradient in the original Wannier90. 
The result is a `OnGrid{<:BrillouinZone}`
with a matrix on each k-point.

```jldoctest orbital_set
julia> G_w = gauge_gradient(U, scheme, brillouin_zone, W90BranchCut);
julia> G_w[brillouin_zone[0, 0, 1]]
4×4 Matrix{ComplexF64}:
       0.0-0.030634im     -1.13757+0.493027im   0.145804-0.00940629im  -0.0627567-0.0302657im
   1.13757+0.493027im          0.0+1.76773im   -0.876742+0.121255im      0.885767+0.370595im
 -0.145804-0.00940629im   0.876742+0.121255im        0.0-0.158586im    -0.0557475+0.0957635im
 0.0627567-0.0302657im   -0.885767+0.370595im  0.0557475+0.0957635im          0.0-0.12397im
```
"""
function gauge_gradient(U::Gauge, scheme::CosScheme, brillouin_zone::B, ::Type{W90BranchCut}) where {B<:BrillouinZone}
    M = gauge_transform(neighbor_basis_integral(scheme), U)
    N = n_band(M)
    c = (n -> center(M, scheme, n, W90BranchCut)).(1:N)

    G = k -> sum(zip(weights(scheme), shells(scheme))) do (w, shell)
        sum(shell) do b
            A = M[k, k+b]
            q = [angle(A[n, n]) + cartesian(b)' * c[n] for n = 1:N]
            R = hcat([A[:, n] * A[n, n]' for n = 1:N]...)
            T = hcat([(A[:, n] / A[n, n]) * q[n] for n = 1:N]...)
            4w * ((R - R') / 2 - (T + T') / 2im) / length(brillouin_zone)
        end
    end

    return map(G, brillouin_zone)
end

"""
    gauge_gradient(U, scheme, brillouin_zone, TruncatedConvolution) 

Gauge gradient for the truncated convolution. 
The result is a `OnGrid{<:BrillouinZone}`
with a matrix on each k-point.

```jldoctest orbital_set
julia> G_t = gauge_gradient(U, scheme, brillouin_zone, TruncatedConvolution);
julia> G_t[brillouin_zone[0, 0, 1]]
4×4 Matrix{ComplexF64}:
        0.0-0.0388009im  -0.00321511+0.158118im    0.160808-0.0451399im  -0.0646812+0.0334942im
 0.00321511+0.158118im           0.0-0.179264im   -0.158035-0.0368048im    0.205487+0.0206957im
  -0.160808-0.0451399im     0.158035-0.0368048im        0.0-0.108496im    -0.070136+0.157672im
  0.0646812+0.0334942im    -0.205487+0.0206957im   0.070136+0.157672im          0.0-0.0640979im
```
"""
function gauge_gradient(U::Gauge, scheme::CosScheme, brillouin_zone::B, ::Type{TruncatedConvolution}) where {B<:BrillouinZone}
    M = gauge_transform(neighbor_basis_integral(scheme), U)
    N = n_band(M)

    ρ̃ = Dict(b => (n -> 1 / N * sum(k -> M[k, k+b][n, n], brillouin_zone)).(1:N) for b in vcat(shells(scheme)...))

    G = k -> sum(zip(weights(scheme), shells(scheme))) do (w, shell)
        sum(shell) do b
            R = hcat([M[k, k+b][:, n] * ρ̃[b][n]' / abs(ρ̃[b][n]) for n = 1:N]...)
            2w * (R - R') / length(brillouin_zone)
        end
    end

    return map(G, brillouin_zone)
end


function all_spread(U, scheme, ::Type{Branch}) where {Branch}
    M = gauge_transform(neighbor_basis_integral(scheme), U)
    N = n_band(M)
    return map(1:N) do n
        spread(M, scheme, n, Branch)
    end
end

function all_center(U, scheme, ::Type{Branch}) where {Branch}
    M = gauge_transform(neighbor_basis_integral(scheme), U)
    N = n_band(M)
    return map(1:N) do n
        center(M, scheme, n, Branch)
    end
end


"""
    ILAOptimizer(scheme)

A vanilla gradient descent optimizer with line search. There are probably faster algorithms,
but this simple optimizer is already out of the scope of this package. 

Example:

```jldoctest orbital_set
julia> optimizer = ILAOptimizer(scheme);
julia> U_optimal = optimizer(U, TruncatedConvolution, FletcherReeves);
julia> M_optimal = gauge_transform(neighbor_basis_integral(scheme), U_optimal);
julia> sum(i->spread(M_optimal, scheme, i, TruncatedConvolution), 1:4)
24.206845069491276
```

Optimization with `W90BranchCut` is not included because it's too much work for 
`Documenter.jl`.
"""
struct ILAOptimizer
    scheme::CosScheme
    meta::Dict{Symbol, Any}
    logger::Dict{Symbol,Vector}
end

ILAOptimizer(scheme::CosScheme) = ILAOptimizer(scheme,
    Dict{Symbol, Any}(),
    Dict{Symbol, Vector}(
        :w90_branch_cut_center => Vector{Matrix{Float64}}(),
        :truncated_convolution_center => Vector{Matrix{Float64}}(),
        :convolutional_center => Vector{Matrix{Float64}}(),
        :w90_branch_cut_spread => Vector{Vector{Float64}}(),
        :truncated_convolution_spread => Vector{Vector{Float64}}(),
        :convolutional_spread => Vector{Vector{Float64}}(),
        :step_size => Vector{Vector{Float64}}()
    ))

scheme(optimizer::ILAOptimizer) = optimizer.scheme
logger(optimizer::ILAOptimizer) = optimizer.logger

"""
    make_step(U, ΔW, α)

Make a step from `log(U)` in the direction of `ΔW` by size `α`.  Keep in mind
that the change in `U` is an approximate since `α ΔW` does not commute with
`log(U)`, and the exponentiation cannot exactly be separated.
"""
function make_step(U::Gauge, ΔW, α)
    brillouin_zone = grid(U)
    new_elements = elements(map(k -> U[k] * cis(-Hermitian(1im * α * ΔW[k])), brillouin_zone))
    V = Gauge{typeof(brillouin_zone)}(brillouin_zone, new_elements)
    return V
end

"""
    line_search(U, f, ∇f, ∇f², α, α_0 = 2)

Perform a line search in the direction of the gradient.
"""
function line_search(U, f, ∇f, ∇f², α; α_0 = 2)
    Ω, ∇Ω = f(U), ∇f(U)
    ∇Ω² = ∇f²(∇Ω)

    while true
        Q = Ω - 0.5α * ∇Ω²
        V = make_step(U, ∇Ω, α)
        Ω_V = f(V)

        discontinuous = ∇Ω² > 1e-7 && α < 1e-3
        discontinuous && return let α = α_0, V = make_step(U, ∇Ω, α), Ω_V = f(V)
            print("✽")
            V, Ω_V, α, ∇Ω²
        end

        Ω_V > Q || return V, Ω_V, α, ∇Ω²
        α /= 2
    end
end

"""
The (tampered) Fletcher Reeves conjugate gradient algorithm.

For an explanation of this conjugate gradient algorithm,
see: http://www.mymathlib.com/optimization/nonlinear/unconstrained/fletcher_reeves.html

The algorithm is adapted alightly since our objective function is not convex.
Line searches can be necessary when the quadratic fit turns out concave.
"""
struct FletcherReeves end

"""
The accelerated gradient descent algortihm. The acceleration starts when the
gradient becomes small enough.
"""
struct AcceleratedGradientDescent end

function (optimizer::ILAOptimizer)(U::Gauge, ::Type{Branch}, ::Type{AcceleratedGradientDescent}; n_iteration = Inf, α_0 = 2, ϵ = 1e-7, logging = false) where {Branch}
    N = optimizer |> scheme |> neighbor_basis_integral |> n_band
    brillouin_zone = grid(U)
    f = U -> sum(all_spread(U, scheme(optimizer), Branch))
    ∇f = U -> Gauge{typeof(brillouin_zone)}(brillouin_zone, elements(gauge_gradient(U, scheme(optimizer), brillouin_zone, Branch)))
    ∇f² = ∇f -> sum(k -> norm(reshape(∇f[k], N^2))^2, brillouin_zone)
    α = α_0
    t = 1
    current_iteration = -1
    while n_iteration > current_iteration
        current_iteration += 1
        X, Ω, α, ∇Ω² = line_search(U, f, ∇f, ∇f², α, α_0 = α_0)
        ∇Ω² / (length(brillouin_zone) * N^2) < ϵ && break
        logging && log(optimizer, U, current_iteration, Ω, ∇Ω², α)
        if ∇Ω² / (length(brillouin_zone) * N^2) > 1e-3
            α = 2α
            U = X
            continue
        end
        t_next = (1 + sqrt(1 + 4t^2)) / 2
        U = X + (t - 1) / t_next * (X - U)
        t = t_next
    end
    return U
end


function (optimizer::ILAOptimizer)(U::Gauge, ::Type{Branch}, ::Type{FletcherReeves}; ϵ = 1e-7, logging = false) where {Branch}
    N = optimizer |> scheme |> neighbor_basis_integral |> n_band
    brillouin_zone = grid(U)
    f = U -> sum(all_spread(U, scheme(optimizer), Branch))
    ∇f = U -> Gauge{typeof(brillouin_zone)}(brillouin_zone, elements(gauge_gradient(U, scheme(optimizer), brillouin_zone, Branch)))
    ∇f² = ∇f -> sum(k -> norm(reshape(∇f[k], N^2))^2, brillouin_zone)
    current_iteration = 0
    h_old, g_old = f(U), ∇f(U)
    g²_old = ∇f²(g_old)
    v_old = Gauge(brillouin_zone, N)
    λ_old = 2

    while true
        g²_old / (length(brillouin_zone) * N^2) < ϵ && break
        g = ∇f(U)
        g² = ∇f²(g)
        α = rem(current_iteration, N^2) == 0 ? 0 : g² / g²_old
        v = Gauge{typeof(brillouin_zone)}(brillouin_zone, elements(map(k -> g[k] + α * v_old[k], brillouin_zone)))
        λ, h = quadratic_fit_1d(λ -> make_step(U, v, λ) |> f)
        if λ > 0 && h < h_old
            U = make_step(U, v, λ)
            h = f(U)
        else
            U, h, λ, _ = line_search(U, f, ∇f, ∇f², 2λ_old)
        end

        logging && log(optimizer, U, current_iteration, h, g², λ)
        λ_old, h_old, v_old, g_old, g²_old = λ, h, v, g, g²
        current_iteration += 1
    end
    return U
end


function log(optimizer, U, current_iteration, Ω, ∇Ω², α)
    println("Iteration: $(current_iteration)")
    println("Ω: $(Ω)")
    println("∇Ω²: $(∇Ω²)")
    println("α: $(α)")
    println()

    N = optimizer |> scheme |> neighbor_basis_integral |> n_band
    append!(logger(optimizer)[:truncated_convolution_spread], [all_spread(U, scheme(optimizer), TruncatedConvolution)])
    append!(logger(optimizer)[:w90_branch_cut_spread], [all_spread(U, scheme(optimizer), W90BranchCut)])
    append!(logger(optimizer)[:truncated_convolution_center], [hcat(all_center(U, scheme(optimizer), TruncatedConvolution)...)])
    append!(logger(optimizer)[:w90_branch_cut_center], [hcat(all_center(U, scheme(optimizer), W90BranchCut)...)])


    # haskey(optimizer.meta, :truncated_convolution_spread)  || (optimizer.meta[:truncated_convolution_spread] = Vector{Vector{Float64}}())
    # haskey(optimizer.meta, :w90_branch_cut_spread)  || (optimizer.meta[:w90_branch_cut_spread] = Vector{Vector{Float64}}())
    # haskey(optimizer.meta, :convolutional_spread)  || (optimizer.meta[:convolutional_spread] = Vector{Vector{Float64}}())

    # haskey(optimizer.meta, :truncated_convolution_center)  || (optimizer.meta[:truncated_convolution_center] = Vector{Matrix{Float64}}())
    # haskey(optimizer.meta, :w90_branch_cut_center)  || (optimizer.meta[:w90_branch_cut_center] = Vector{Matrix{Float64}}())
    # haskey(optimizer.meta, :convolutional_center)  || (optimizer.meta[:convolutional_center] = Vector{Matrix{Float64}}())


    mod(current_iteration, 20) == 0 || return
    haskey(optimizer.meta, :ũ) || return
    ũ = set_gauge(optimizer.meta[:ũ], U)
    r̃2 = optimizer.meta[:r̃2]
    ρ̃ = reciprocal_densities(ũ)

    convolutional_center = zeros(3, length(ρ̃))
    convolutional_spread = zeros(length(ρ̃))
    for i in 1:length(ρ̃)
        c, σ² = center_spread(ρ̃[i], r̃2)
        convolutional_center[:, i] = c
        convolutional_spread[i] = σ²
    end
    append!(logger(optimizer)[:convolutional_center], [convolutional_center])
    append!(logger(optimizer)[:convolutional_spread], [convolutional_spread])
end

