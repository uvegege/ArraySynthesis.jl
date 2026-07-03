abstract type AbstractDirection end

struct ThetaDirection{T} <: AbstractDirection
    θ::T
end

struct UVDirection{T} <: AbstractDirection
    u::T
    v::T
end

const Direction{T} = Union{ThetaDirection{T}, UVDirection{T}}

struct Region{P}
    points::Vector{P}
    name::Symbol
end

abstract type AbstractPatternItem{P,T} end

struct Beam{P, T} <: AbstractPatternItem{P,T}
    direction::P
    gain::T
end

struct ShapedBeam{P, T} <: AbstractPatternItem{P,T}
    region::Region{P}
    target::Vector{T}
    ripple::T
end

struct NullPoint{P, T} <: AbstractPatternItem{P,T}
    direction::P
    level::T
end

struct Nulls{P,T} <: AbstractPatternItem{P,T}
    points::Vector{NullPoint{P,T}}
end

struct SideLobeRegion{P, T} <: AbstractPatternItem{P,T}
    region::Region{P}
    upper::Vector{T}
end

struct Pattern{P, T}
    beams::Vector{Beam{P, T}}
    shaped_beams::Vector{ShapedBeam{P, T}}
    null_directions::Vector{NullPoint{P, T}}
    sidelobe_regions::Vector{SideLobeRegion{P, T}}
end

#=
additems!(p) = p

function additems!(p, item, rest...)
    additem!(p, item)
    additems!(p, rest...)
end

additem!(p, x::Beam) = push!(p.beams, x)
additem!(p, x::ShapedBeam) = push!(p.shaped_beams, x)
additem!(p, x::NullPoint) = push!(p.null_directions, x)
additem!(p, x::Nulls) = append!(p.null_directions, x.points)
additem!(p, x::SideLobeRegion) = push!(p.sidelobe_regions, x)

function pattern(items::AbstractPatternItem{P,T}...) where {P,T}
    p = Pattern{P,T}(
        Beam{P,T}[],
        ShapedBeam{P,T}[],
        NullPoint{P,T}[],
        SideLobeRegion{P,T}[]
    )

    additems!(p, items...)

    return p
end
=#

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
            error("¿?")
        end
    end
    return Pattern(beams, shaped_beams, null_points, sidelobe_regions)
end

struct ClosedInterval{T}
    a::T
    b::T
end

..(a, b) = ClosedInterval(a, b)

const ° = pi / 180
struct dB end
Base.:*(x, ::Type{dB}) = 10^(x/20)

direction(x::ThetaDirection) = x
direction(x::UVDirection) = x
direction(x::Number) = ThetaDirection(x)
direction(x::NamedTuple{(:u, :v)}) = UVDirection(x.u, x.v)

θ(x) = ThetaDirection(x)
uv(θ, ϕ) = UVDirection(sin(θ)*cos(ϕ), sin(θ)*sin(ϕ))
@inline w(u, v) = sqrt(1 - u^2 - v^2)

beam(dir; gain = 1.0) = Beam(direction(dir), gain)
null(dir; level = -60.0dB) = NullPoint(direction(dir), level)
nulls(dirs; level = -60.0dB) = Nulls([null(d; level) for d in dirs])


function region(r::ClosedInterval; npoints = 60, name = :region)
    points = [θ(x) for x in range(r.a, r.b, npoints)]
    return Region(points, name)
end

function region(r::ClosedInterval, step; name = :region)
    points = [θ(x) for x in r.a:step:r.b]
    return Region(points, name)
end


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

function join_regions(x, y)
    Region(vcat(x.points, y.points), x.name)
end
