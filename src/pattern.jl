"""
Abstract supertype for all direction parameterizations.
"""
abstract type AbstractDirection end

"""
    ThetaDirection(θ)

One-dimensional angular direction. Angles are in radians; use `°` for degrees.
"""
struct ThetaDirection{T} <: AbstractDirection
    θ::T
end

"""
    UVDirection(u, v)

Direction in the visible `(u, v)` plane.
"""
struct UVDirection{T} <: AbstractDirection
    u::T
    v::T
end

"""
Union alias for `ThetaDirection` and `UVDirection`.
"""
const Direction{T} = Union{ThetaDirection{T}, UVDirection{T}}

"""
    Region(points, name)

Named collection of sampled directions used for masks, shaped beams, or
integration regions.
"""
struct Region{P}
    points::Vector{P}
    name::Symbol
end

abstract type AbstractPatternItem{P,T} end

"""
    Beam(direction, gain)

Main-beam constraint at `direction` with the requested linear `gain`.
Prefer the `beam` constructor in user code.
"""
struct Beam{P, T} <: AbstractPatternItem{P,T}
    direction::P
    gain::T
end

"""
    ShapedBeam(region, target, ripple)

Shaped-beam mask over a sampled `Region`. Prefer `shaped_beam` in user code.
"""
struct ShapedBeam{P, T} <: AbstractPatternItem{P,T}
    region::Region{P}
    target::Vector{T}
    ripple::T
end

"""
    NullPoint(direction, level)

Null constraint at one direction. Prefer `null` in user code.
"""
struct NullPoint{P, T} <: AbstractPatternItem{P,T}
    direction::P
    level::T
end

"""
    Nulls(points)

Collection of null constraints. Prefer `nulls` in user code.
"""
struct Nulls{P,T} <: AbstractPatternItem{P,T}
    points::Vector{NullPoint{P,T}}
end

"""
    SideLobeRegion(region, upper)

Sidelobe upper mask over a sampled `Region`. Prefer `sidelobes` in user code.
"""
struct SideLobeRegion{P, T} <: AbstractPatternItem{P,T}
    region::Region{P}
    upper::Vector{T}
end

"""
    Pattern(beams, shaped_beams, null_directions, sidelobe_regions)

Container for all pattern constraints passed to `synthesize`.
Prefer the flexible `pattern` constructor in user code.
"""
struct Pattern{P, T}
    beams::Vector{Beam{P, T}}
    shaped_beams::Vector{ShapedBeam{P, T}}
    null_directions::Vector{NullPoint{P, T}}
    sidelobe_regions::Vector{SideLobeRegion{P, T}}
end

"""
    pattern(items...)

Build a `Pattern` from any mixture of `beam`, `shaped_beam`, `null`, `nulls`,
and `sidelobes` items.

This is the main user-facing pattern builder: items can be written in the order
that is most natural for the synthesis problem, and `pattern` sorts them into
the internal beam, shaped-beam, null, and sidelobe collections. All items must use
compatible direction and numeric types.

# Examples

```julia
p = pattern(
    beam(0°),
    nulls([-30°, 30°]),
    sidelobes(region(-90°..-10°, 1°), -25dB),
    sidelobes(region(10°..90°, 1°), theta_ramp(10°, -30dB, 90°, -20dB)),
)
```
"""
function pattern(items::AbstractPatternItem{P, T}...) where {P, T}
    beams = Beam{P, T}[]
    shaped_beams = ShapedBeam{P, T}[]
    null_points = NullPoint{P, T}[]
    sidelobe_regions = SideLobeRegion{P, T}[]
    for item in items
        if item isa Beam{P, T}
            push!(beams, item)
        elseif item isa ShapedBeam{P, T}
            push!(shaped_beams, item)
        elseif item isa Nulls{P, T}
            append!(null_points, item.points)
        elseif item isa NullPoint{P, T}
            push!(null_points, item)
        elseif item isa SideLobeRegion{P, T}
            push!(sidelobe_regions, item)
        else
            error("wrong item")
        end
    end
    return Pattern(beams, shaped_beams, null_points, sidelobe_regions)
end

"""
    ClosedInterval(a, b)
    a..b

Closed interval helper used to create angular regions.
"""
struct ClosedInterval{T}
    a::T
    b::T
end

"""
    a..b

Create a `ClosedInterval(a, b)`.
"""
..(a, b) = ClosedInterval(a, b)

const ° = pi / 180
struct dB end
Base.:*(x, ::Type{dB}) = 10^(x/20)

"""
    direction(x)

Convert numbers and `(u = ..., v = ...)` named tuples into direction objects.
"""
direction(x::ThetaDirection) = x
direction(x::UVDirection) = x
direction(x::Number) = ThetaDirection(x)
direction(x::NamedTuple{(:u, :v)}) = UVDirection(x.u, x.v)

"""
    θ(x)

Construct a `ThetaDirection`.
"""
θ(x) = ThetaDirection(x)

"""
    uv(θ, ϕ)

Convert spherical angles to a `UVDirection`.
"""
uv(θ, ϕ) = UVDirection(sin(θ)*cos(ϕ), sin(θ)*sin(ϕ))
@inline w(u, v) = sqrt(1 - u^2 - v^2)

"""
    beam(dir; gain = 1.0)

Create a main-beam constraint at `dir`.
"""
beam(dir; gain = 1.0) = Beam(direction(dir), gain)

"""
    null(dir; level = -60.0dB)

Create a null constraint at `dir` with the given maximum level.
"""
null(dir; level = -60.0dB) = NullPoint(direction(dir), level)

"""
    nulls(dirs; level = -60.0dB)

Create null constraints for all directions in `dirs`.
"""
nulls(dirs; level = -60.0dB) = Nulls([null(d; level) for d in dirs])


"""
    region(interval; npoints = 60, name = :region)
    region(interval, step; name = :region)

Sample a one-dimensional angular interval into a `Region` of `ThetaDirection`
points.
"""
function region(r::ClosedInterval; npoints = 60, name = :region)
    points = [θ(x) for x in range(r.a, r.b, npoints)]
    return Region(points, name)
end

function region(r::ClosedInterval, step; name = :region)
    points = [θ(x) for x in r.a:step:r.b]
    return Region(points, name)
end


"""
    shaped_beam(region, target; ripple = 1.0dB, normalize = false)

Create a shaped-beam constraint over `region`.

`target` may be a scalar, a vector with one value per region point, or a
function evaluated at each point. When `normalize` is true, target values are
scaled by their maximum absolute value.
"""
function shaped_beam(region::Region, target::F; ripple = 1.0dB, normalize = true) where F <: Function
    values = target.(region.points)
    if normalize
        m = maximum(abs, values)
        m == 0 && error("error in shaped normalization")
        values ./= m
    end
    return ShapedBeam(region, values, ripple)
end

function shaped_beam(region::Region, target::AbstractVector; ripple = 1.0dB, normalize = false)
    length(target) == length(region.points) || error("Target length must match region points.")
    values = collect(target)
    if normalize
        m = maximum(abs, values)
        m == 0 && error("error in shaped normalization")
        values ./= m
    end
    return ShapedBeam(region, values, ripple)
end

function shaped_beam(region::Region, target::Number; ripple = 1.0dB)
    values = fill(target, length(region.points))
    return ShapedBeam(region, values, ripple)
end


"""
    sidelobes(region, upper = -20dB)

Create an upper sidelobe mask over `region`.

`upper` may be a scalar, a vector with one value per region point, or a
function evaluated at each point.
"""
function sidelobes(region::Region, upper = -20dB)
    values = evaluate_mask(upper, region)
    return SideLobeRegion(region, values)
end

sidelobes(interval::ClosedInterval, u) = sidelobes(region(interval; npoints = 60), u)


evaluate_mask(x::Number, region) = fill(x, length(region.points))
evaluate_mask(f::F, region) where F <: Function = [f(p) for p in region.points]
function evaluate_mask(v::AbstractVector, region)
    length(v) == length(region.points) || error("Mask length must match region points")
    return collect(v)
end

"""
    outside(intervals...; limits = -90°..90°)

Return the complementary closed intervals inside `limits`.
"""
function outside(r::ClosedInterval; limits = -90°..90°)
    regions = ClosedInterval[]
    if r.a > limits.a
        push!(regions, limits.a..r.a)
    end
    if r.b < limits.b
        push!(regions, r.b..limits.b)
    end
    return regions
end

outside(intervals::ClosedInterval...; limits = -90°..90°) = outside(ClosedInterval[intervals...]; limits = limits)

function outside(intervals::AbstractVector{<:ClosedInterval}; limits = -90°..90°)
    isempty(intervals) && return [limits]

    sorted = sort(intervals; by = r -> r.a)

    for r in sorted
        r.a < limits.a && error("Interval $(r) starts below limits $(limits).")
        r.b > limits.b && error("Interval $(r) ends above limits $(limits).")
        r.a > r.b && error("Invalid interval $(r).")
    end

    for i in 1:length(sorted)-1
        r1, r2 = sorted[i], sorted[i+1]
        r1.b > r2.a && error("Intervals $(r1) and $(r2) overlap.")
    end

    regions = ClosedInterval[]
    if first(sorted).a > limits.a
        push!(regions, limits.a..first(sorted).a)
    end

    for i in 1:length(sorted)-1
        r1, r2 = sorted[i], sorted[i+1]
        if r1.b < r2.a
            push!(regions, r1.b..r2.a)
        end
    end

    if last(sorted).b < limits.b
        push!(regions, last(sorted).b..limits.b)
    end

    return regions
end

"""
    join_regions(x, y)

Concatenate two sampled regions into a single `Region`.
"""
function join_regions(x, y)
    Region(vcat(x.points, y.points), x.name)
end
