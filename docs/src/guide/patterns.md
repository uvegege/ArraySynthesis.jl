# Pattern and Regions

`pattern` is the main user-facing specification layer in ArraySynthesis. It collects
beams, shaped beams, nulls, and sidelobe masks in whatever order is most natural for
the problem:

```julia
p = pattern(
    beam(0°),
    nulls([-30°, 30°]),
    sidelobes(region(-90°..-10°, 1°), -25dB),
    sidelobes(region(10°..90°, 1°), theta_ramp(10°, -30dB, 90°, -20dB)),
)
```

The individual items are sorted internally into a `Pattern`, which is then passed to
[`synthesize`](@ref).

## Direction Systems

One-dimensional cuts use `ThetaDirection`:

```julia
θ0 = θ(15°)
θ0 = ThetaDirection(15°)
```

Planar and conformal-array masks usually use direction cosines:

```julia
p = UVDirection(0.2, -0.1)
p = direction((u = 0.2, v = -0.1))
```

The helper `uv(θ, ϕ)` converts spherical angles into a `UVDirection`.

## One-Dimensional Regions

Use closed intervals and `region` for angular cuts:

```julia
main = region(-10°..10°, 0.5°)
sll_left = region(-90°..-15°, 1°)
sll_right = region(15°..90°, 1°)
```

`outside` builds complementary intervals, which is convenient for linear arrays:

```julia
main_lobe = -12°..12°
sll = join_regions(region.(outside(main_lobe), 1°)...)
```

## Pattern Items

Use `beam` for exact gain constraints:

```julia
beam(0°)
beam((u = 0.1, v = 0.0), gain = 1.0)
```

Use `shaped_beam` when the main region has a mask rather than a single direction:

```julia
flat = shaped_beam(region(-10°..10°, 1°), 1.0, ripple = -0.5dB)
tilted = shaped_beam(region(-10°..10°, 1°), theta_ramp(-10°, 0.8, 10°, 1.0))
```

Use `sidelobes` for upper masks:

```julia
sidelobes(region(15°..90°, 1°), -25dB)
sidelobes(region(15°..90°, 1°), theta_ramp(15°, -35dB, 90°, -20dB))
```

Use `null` and `nulls` for deep constraints at isolated directions:

```julia
null(30°)
nulls([-30°, 30°], level = -60dB)
```

## Planar Complementary Regions

For planar and conformal arrays, the easiest specification is often geometric:
draw the beam or protected footprint in the `(u, v)` plane, then synthesize over
the complementary visible region. That is the role of [`visible_region`](@ref).

```julia
beam_shape = Ellipse(0.25, 0.15, (0.0, 0.0))
beam_region = region(beam_shape; step = 0.02)
sll_region = visible_region(beam_shape; step = 0.02, bandpass = 0.05,
                             filtered = true)

p = pattern(
    shaped_beam(beam_region, 1.0, ripple = -0.5dB),
    sidelobes(sll_region, -25dB),
)
```

`visible_region(shapes...)` samples the `(u, v)` visible plane and removes all
directions inside the supplied shapes. `bandpass` expands the removed zones to
leave a guard band between the shaped/protected region and the sidelobe mask.

Available shapes include `Circle`, `Ellipse`, `Rectangle`, `Polygon`, `Moonlike`,
and convenience constructors such as `rhombus` and `triangle`.
