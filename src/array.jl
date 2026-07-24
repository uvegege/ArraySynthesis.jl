"""
    ArrayGeometry(positions, dim)

Array geometry with one column per element in the `3 x N` `positions` matrix.
Coordinates are expressed in wavelengths, and `dim` is the effective support
dimension of the array.
"""
struct ArrayGeometry{T}
    positions::Matrix{T}
    dim::Int
end

"""
    SymmetricArray(positions, dim)

Array represented by only one element from each conjugate-symmetric pair.
Use `materialize` to expand it into the corresponding full `ArrayGeometry`.
"""
struct SymmetricArray{T}
    positions::Matrix{T}
    dim::Int
end

function support_dim(p::AbstractMatrix)
    return count(i -> any(!iszero, p[i, :]), 1:3)
end

"""
    is_origin(p)

Return whether point `p` is the array origin.
"""
function is_origin(p)
    return all(iszero, p)
end

"""
    is_representative(p)

Return whether point `p` is kept as the representative of a symmetric pair.
"""
function is_representative(p)
    is_origin(p) && return true
    i = findfirst(!iszero, p)
    return p[i] > 0
end

"""
    uniform_linear_array(N; d = 0.5, axis = :x, centered = true)

Create a uniform linear array with `N` elements spaced by `d` wavelengths.
"""
function uniform_linear_array(N::Integer; d=0.5, axis=:x, centered=true)
    p = zeros(Float64, 3, N)
    x = centered ? (collect(0:N-1) .- (N-1)/2) .* d : collect(0:N-1) .* d
    dim = axis === :x ? 1 : axis === :y ? 2 : 3
    p[dim, :] .= x
    return ArrayGeometry(p, 1)
end

"""
    linear_array(x; axis = :x)

Create a linear array with element coordinates given by `x`.
"""
function linear_array(x::AbstractVector; axis=:x)
    p = zeros(eltype(x), 3, length(x))
    dim = axis === :x ? 1 : axis === :y ? 2 : 3
    p[dim, :] .= x
    return ArrayGeometry(p, 1)
end

"""
    symmetric_linear_array(N; d = 0.5, axis = :x)

Create a conjugate-symmetric representation of a centered uniform linear array.
"""
function symmetric_linear_array(N::Integer; d=0.5, axis=:x)
    full = uniform_linear_array(N; d, axis, centered=true).positions
    p, dim = symmetrize(full)
    return SymmetricArray(p, dim)
end

"""
    planar_array(Nx, Ny; dx = 0.5, dy = 0.5, centered = true)

Create a rectangular planar array with `Nx * Ny` elements.
"""
function planar_array(Nx::Integer, Ny::Integer; dx=0.5, dy=0.5, centered=true)
    xs = collect(0:Nx-1) .* dx
    ys = collect(0:Ny-1) .* dy
    centered && (xs .-= sum(xs) / length(xs); ys .-= sum(ys) / length(ys))
    p = zeros(Float64, 3, Nx*Ny)

    k = 1
    for y in ys
        for x in xs
            p[:, k] .= (x, y, 0.0)
            k += 1
        end
    end

    return ArrayGeometry(p, 2)
end

"""
    triangular_array(Nx, Ny; d = 0.5, centered = true)

Create a planar triangular-lattice array.
"""
function triangular_array(Nx::Integer, Ny::Integer; d=0.5, centered=true)
    p = zeros(Float64, 3, Nx*Ny)

    k = 1
    for j in 0:Ny-1
        for i in 0:Nx-1
            p[1, k] = d * (i + 0.5 * isodd(j))
            p[2, k] = d * sqrt(3) / 2 * j
            k += 1
        end
    end

    centered && (p .-= sum(p; dims=2) ./ size(p, 2))

    return ArrayGeometry(p, 2)
end

"""
    hexagonal_array(order; d = 0.5)

Create a centered hexagonal array with the given lattice `order`.
"""
function hexagonal_array(order::Integer; d=0.5)
    N = 1 + 3 * order * (order + 1)
    p = zeros(Float64, 3, N)

    k = 1
    for q in -order:order
        for r in max(-order, -q-order):min(order, -q+order)
            p[1, k] = d * (q + r/2)
            p[2, k] = d * sqrt(3) / 2 * r
            k += 1
        end
    end

    return ArrayGeometry(p, 2)
end

"""
    circular_array(radii, elements)

Create concentric circular rings, with `elements[i]` elements on `radii[i]`.
"""
function circular_array(radii::AbstractVector, elements::AbstractVector)
    length(radii) == length(elements) || throw(ArgumentError("radii and elements must have the same length"))

    N = sum(elements)
    p = zeros(Float64, 3, N)

    k = 1
    for (ρ, M) in zip(radii, elements)
        if ρ == 0
            p[:, k] .= (0.0, 0.0, 0.0)
            k += 1
        else
            for n in 0:M-1
                φ = 2π * n / M
                p[1, k] = ρ * cos(φ)
                p[2, k] = ρ * sin(φ)
                k += 1
            end
        end
    end

    return ArrayGeometry(p, 2)
end

"""
    symmetric_planar_array(Nx, Ny; dx = 0.5, dy = 0.5, centered = true)

Create a conjugate-symmetric representative of a rectangular planar array.
"""
function symmetric_planar_array(Nx::Integer, Ny::Integer; dx=0.5, dy=0.5, centered=true)
    full = planar_array(Nx, Ny; dx, dy, centered).positions
    p, dim = symmetrize(full)
    return SymmetricArray(p, dim)
end

"""
    symmetric_triangular_array(Nx, Ny; d = 0.5, centered = true)

Create a conjugate-symmetric representative of a triangular-lattice array.
"""
function symmetric_triangular_array(Nx::Integer, Ny::Integer; d=0.5, centered=true)
    full = triangular_array(Nx, Ny; d, centered).positions
    sym_positions, dim = symmetrize(full)
    return SymmetricArray(sym_positions, dim)
end

"""
    symmetric_hexagonal_array(order; d = 0.5)

Create a conjugate-symmetric representative of a hexagonal array.
"""
function symmetric_hexagonal_array(order::Integer; d=0.5)
    full = hexagonal_array(order; d).positions
    sym_positions, dim = symmetrize(full)
    return SymmetricArray(sym_positions, dim)
end

"""
    symmetric_circular_array(radii, elements)

Create a conjugate-symmetric representative of concentric circular rings.
"""
function symmetric_circular_array(radii::AbstractVector, elements::AbstractVector)
    full = circular_array(radii, elements).positions
    sym_positions, dim = symmetrize(full)
    return SymmetricArray(sym_positions, dim)
end

"""
    symmetrize(full)

Keep one representative per conjugate-symmetric pair in a full `3 x N`
position matrix and return `(positions, dim)`.
"""
function symmetrize(full)
    dim = support_dim(full)
    keep = filter(is_representative, eachcol(full))
    p = zeros(eltype(full), 3, length(keep))
    for k in eachindex(keep)
        p[:, k] .= keep[k]
    end
    return p, dim
end

"""
    materialize(array::SymmetricArray)

Expand a symmetric representative array into its full physical `ArrayGeometry`.
"""
function materialize(a::SymmetricArray)
    p = a.positions
    N = sum(col -> is_origin(col) ? 1 : 2, eachcol(p))
    full = zeros(eltype(p), 3, N)

    k = 1
    for col in eachcol(p)
        full[:, k] .= col
        k += 1
        is_origin(col) && continue
        full[:, k] .= -col
        k += 1
    end
    return ArrayGeometry(full, a.dim)
end
