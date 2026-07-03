using ArraySynthesis
using ArraySynthesis: °, dB
using GLMakie
using LinearAlgebra
using HiGHS


# Resultados_2-style: 16x16 symmetric planar, diamond beam at (14°,12°), null region.
array_full = planar_array(16, 16, dx = 0.5, dy = 0.5)
array = SymmetricArray(symmetrize(array_full.positions)...)

focus = uv(14°, 12°)
beam_shape = rhombus((focus.u, focus.v), 16°)
beam_region = region(beam_shape, step = 2°)

null_dir = uv(-55°, -55°)
null_shape = Circle(0.3, (null_dir.u, null_dir.v))
null_region = region(null_shape, step = 2°)
sl_region = visible_region(beam_shape, null_shape; step = 4°, bandpass = 0.12)

p = pattern(shaped_beam(beam_region, 1.0, ripple = -1dB),
    sidelobes(null_region, -50dB))
obj = MinSLL([sl_region])

result1 = synthesize(array, p, obj, ConjugateSymmetricWeights(), LP(), HiGHS.Optimizer)
result2 = synthesize(array, p, obj, ProgressivePhaseAmplitude(focus), LP(), HiGHS.Optimizer)

U = V = collect(-1.0:0.02:1.0)
dirs = [UVDirection(u, v) for u in U, v in V]
coef = ConjugateSymmetricWeights()
AF = array_factor(array, coef, result1.weights, dirs)
af_vals = reshape(20 .* log10.(abs.(getindex.(AF, 1)) .+ 1e-12), length(U), length(V))
AF2 = array_factor(array, ProgressivePhaseAmplitude(focus), result2.weights, dirs)
af_vals2 = reshape(20 .* log10.(abs.(getindex.(AF2, 1)) .+ 1e-12), length(U), length(V))

fig = Figure()
ax = Axis(fig[1,1], xlabel="u", ylabel="v")
image!(ax, (-1,1), (-1,1), af_vals, colorrange=(-40, 0), colormap = parula_cm)
Colorbar(fig[1,3], colorrange=(-40, 0)); fig
ax2 = Axis3(fig[1,2])
surface!(ax2, U, V, map(x->max(-50, x), af_vals), colormap = parula_cm)
ax = Axis(fig[2,1], xlabel="u", ylabel="v")
image!(ax, (-1,1), (-1,1), af_vals2, colorrange=(-40, 0), colormap = parula_cm)
Colorbar(fig[2,3], colorrange=(-40, 0)); fig
ax2 = Axis3(fig[2,2])
surface!(ax2, U, V, map(x->max(-50, x), af_vals2), colormap = parula_cm)
fig

