"""
    RegionShape

Abstract type for two-dimensional shapes sampled in the visible `(u, v)` plane.
"""
abstract type RegionShape end

"""Circular region shape in the `(u, v)` plane."""
struct Circle{T} <: RegionShape
    radius::T
    center::Tuple{T, T}
end

"""Elliptic region shape in the `(u, v)` plane."""
struct Ellipse{T} <: RegionShape
    a::T
    b::T
    center::Tuple{T, T}
end

"""Axis-aligned rectangular region shape in the `(u, v)` plane."""
struct Rectangle{T} <: RegionShape
    lx::T
    ly::T
    center::Tuple{T, T}
end

"""Polygonal region shape defined by vertex coordinates `u` and `v`."""
struct Polygon{T} <: RegionShape
    u::Vector{T}
    v::Vector{T}
end

"""Moon-like shape defined as one disk with a second disk removed."""
struct Moonlike{T} <: RegionShape
    r1::T
    r2::T
    center1::Tuple{T, T}
    center2::Tuple{T, T}
end

"""
    rhombus(center, d)

Create a rhombus-shaped `Polygon` centered at `center`.
"""
function rhombus(center, d)
    u0, v0 = center
    return Polygon([u0, u0 + d, u0, u0 - d],[v0 + d, v0, v0 - d, v0])
end

"""
    triangle(center; base, height)

Create a triangular `Polygon` centered at `center`.
"""
function triangle(center; base, height)
    u0, v0 = center
    return Polygon([u0 - base/2, u0 + base/2, u0], [v0 - height/2, v0 - height/2, v0 + height/2])
end

bbox(c::Circle) = (c.center[1] - c.radius, c.center[1] + c.radius, c.center[2] - c.radius, c.center[2] + c.radius)
bbox(e::Ellipse) = (e.center[1] - e.a, e.center[1] + e.a, e.center[2] - e.b, e.center[2] + e.b)
bbox(r::Rectangle) = (r.center[1] - r.lx,r.center[1] + r.lx,r.center[2] - r.ly,r.center[2] + r.ly)
bbox(p::Polygon) = (minimum(p.u),maximum(p.u),minimum(p.v),maximum(p.v))

function bbox(m::Moonlike)
    u1min, u1max, v1min, v1max = bbox(Circle(m.r1, m.center1))
    u2min, u2max, v2min, v2max = bbox(Circle(m.r2, m.center2))
    return (min(u1min, u2min),max(u1max, u2max),min(v1min, v2min),max(v1max, v2max))
end

function distance(shape::Circle, u::Real, v::Real)
    uo, vo = shape.center
    return sqrt((u - uo)^2 + (v - vo)^2) - shape.radius
end

function distance(shape::Ellipse, u::Real, v::Real)
    uo, vo = shape.center
    return sqrt(((u - uo) / shape.a)^2 + ((v - vo) / shape.b)^2) - 1.0
end

function distance(shape::Rectangle, u::Real, v::Real)
    uo, vo = shape.center

    dx = abs(u - uo) - shape.lx
    dy = abs(v - vo) - shape.ly

    ax = max(dx, 0)
    ay = max(dy, 0)

    outside = sqrt(ax^2 + ay^2)
    inside = min(max(dx, dy), 0)

    return outside + inside
end

function distance(shape::Moonlike, u::Real, v::Real)
    uo, vo = shape.center1
    uc, vc = shape.center2

    d1 = sqrt((u - uo)^2 + (v - vo)^2) - shape.r1
    d2 = shape.r2 - sqrt((u - uc)^2 + (v - vc)^2)

    return max(d1, d2)
end

function distance(shape::Polygon, u::Real, v::Real)
    vx = shape.u
    vy = shape.v

    N = length(vx)

    px = u
    py = v

    dx0 = px - vx[1]
    dy0 = py - vy[1]
    d = dx0^2 + dy0^2

    s = 1.0
    j = N

    for i in 1:N
        vi_x = vx[i]
        vi_y = vy[i]
        vj_x = vx[j]
        vj_y = vy[j]

        ex = vj_x - vi_x
        ey = vj_y - vi_y

        wx = px - vi_x
        wy = py - vi_y

        edot = ex*ex + ey*ey
        t = edot == 0 ? 0.0 : clamp((wx*ex + wy*ey) / edot, 0.0, 1.0)

        bx = wx - ex * t
        by = wy - ey * t

        d = min(d, bx*bx + by*by)
        c1 = (py >= vi_y)
        c2 = (py <  vj_y)
        cross = ex*wy > ey*wx

        if (c1 == c2 == cross)
            s = -s
        end

        j = i
    end

    return s * sqrt(d)
end

isinside(shape::T, u, v) where T = distance(shape, u, v) <= 0

filter_inside(shape, uv::UVDirection) = isinside(shape, uv.u, uv.v)

function sample(shape::RegionShape; step::Real, name::Symbol = :region)
    umin, umax, vmin, vmax = bbox(shape)
    points = (UVDirection(u, v) for (u, v) in Iterators.product(umin:step:umax, vmin:step:vmax))
    filtered_points = Iterators.filter(uv -> (1 - uv.u^2 - uv.v^2) > 0, points)
    filtered_points = Iterators.filter(Base.Fix1(filter_inside, shape), filtered_points)
    return Region(vec(collect(filtered_points)), name)
end

"""
    region(shape::RegionShape; step, name = :region)

Sample a two-dimensional shape into a `Region` of `UVDirection` points.
"""
region(shape::RegionShape; step::Real, name::Symbol = :region) = sample(shape; step, name = name)

"""
    visible_region(shapes...; step = 1°, bandpass = 0.05, filtered = false)

Sample the complementary visible-region mask in the `(u, v)` plane.

This is the main helper for planar and conformal-array specifications where
the desired beam or protected zones are easier to draw as shapes. Each supplied
shape is excluded from the sampled visible plane, so the returned `Region`
contains the directions outside all supplied shapes. Use it to define sidelobe
or background regions complementary to shaped beams, null-protection zones, or
other user-defined footprints.

`bandpass` expands the excluded zone around each shape, creating a guard band
between the kept visible region and the removed shapes. When `filtered` is
true, only physically visible points satisfying `u^2 + v^2 < 1` are returned;
when false, the full square `[-1, 1] x [-1, 1]` grid is sampled.

# Example

```julia
main = Ellipse(0.25, 0.15, (0.0, 0.0))
sl = visible_region(main; step = 0.02, bandpass = 0.05, filtered = true)
p = pattern(shaped_beam(region(main; step = 0.02), 1.0),
            sidelobes(sl, -25dB))
```
"""
function visible_region(items::RegionShape...; step = 1°, bandpass = 0.05, filtered = false)
    points = (UVDirection(u, v) for (u, v) in Iterators.product(-1.0:step:1.0, -1.0:step:1.0))
    if filtered
        filtered_points = Iterators.filter(uv -> (1 - uv.u^2 - uv.v^2) > 0, points)
    else
        filtered_points = points
    end
    for item in items
        filtered_points = Iterators.filter(uv -> distance(item, uv.u, uv.v) - bandpass >= 0, filtered_points)
    end
    return Region(vec(collect(filtered_points)), :visible_region)
end
