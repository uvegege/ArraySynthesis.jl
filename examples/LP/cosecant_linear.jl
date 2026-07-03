using ArraySynthesis
using ArraySynthesis: °, dB
using GLMakie
using LinearAlgebra
using HiGHS
using Clarabel


sl_region1 = region(-90°..(-17.5°), 1°)
sl_region2 = region(-17.5°..0°, 1°)
sl_region3 = region(36.87°..90°, 1°)
sb_region = region(5.74°..30°, 1°)

csc_target = map(p -> 0.1 / sin(p.θ), sb_region.points)

p = pattern(shaped_beam(sb_region, csc_target, ripple = -0.6dB, normalize = true),
    sidelobes(sl_region1, -25dB),
    sidelobes(sl_region2, -18dB),
    sidelobes(sl_region3, -25dB))

array = uniform_linear_array(30, d = 0.5)
obj = MinSLL([sl_region1, sl_region2, sl_region3])
coef = ComplexWeights()
result = synthesize(array, p, obj, coef, LP(), HiGHS.Optimizer)

theta_vals = -π/2:0.0005:π/2
af_vals = [abs(array_factor(array, coef, result.weights, [ThetaDirection(θ)])[1]) for θ in theta_vals]
af_db = 20 .* log10.(max.(af_vals, 1e-12))

fig = Figure()
ax = Axis(fig[1,1], xlabel="sin(θ)", ylabel="|AF| (dB)")
lines!(ax, sin.(theta_vals), af_db, linewidth=2)
xlims!(-1, 1); fig
ylims!(-80, 5)
fig