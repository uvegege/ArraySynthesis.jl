using ArraySynthesis
using ArraySynthesis: °, dB
using GLMakie
using LinearAlgebra
using HiGHS
using Mosek, MosekTools

# Semi-axes: a=sin(12°)≈0.208 in u, b=sin(1°)≈0.017 in v (narrow-fan beam)
array = planar_array(13, 13, dx = 0.5, dy = 0.5)

beam_shape = Ellipse(12°, 1°, (0.0, 0.0))
guard_shape = Ellipse(12° + 0.2, 1° + 0.2, (0.0, 0.0))
null_shape = Circle(8°, (0.0, -0.5))
beam_region = region(beam_shape, step = 1°)
sl_region  = visible_region(guard_shape, null_shape; step = 4°, bandpass = 0.0, filtered = false)
null_region = region(null_shape, step = 2°)

p = pattern(shaped_beam(beam_region, 1.0, ripple = -1dB),
    sidelobes(sl_region, -22.4dB),
    sidelobes(null_region, -50dB))

obj = IterativeReweightedL1(max_iter = 15)
coef = ComplexWeights()
result = synthesize(array, p, obj, coef, SOCP(), Mosek.Optimizer)

active = abs.(result.weights) .> 1e-3

begin
fig = Figure()
ax = Axis(fig[1,1], xlabel="x/λ", ylabel="y/λ", title="Active: $(sum(active)) / $(length(active))")
scatter!(ax, array.positions[1, .!active], array.positions[2, .!active], color=:red, marker=:x, markersize=6)
scatter!(ax, array.positions[1,   active], array.positions[2,   active], color=:green, markersize=8)
fig

ax2 = Axis(fig[1,2], xlabel = "u", ylabel = "v", title = "Regions", aspect = DataAspect())
scatter!(ax2, map(uv->(uv.u, uv.v), sl_region.points))
scatter!(ax2, map(uv->(uv.u, uv.v), beam_region.points))
scatter!(ax2, map(uv->(uv.u, uv.v), null_region.points))
U = collect(-1.0:0.01:1.0)
V = collect(-1.0:0.01:1.0)
dirs = [UVDirection(u,v) for u in U, v in V]
AF = array_factor(array, coef, result.weights, dirs)
af_vals = reshape(20 .* log10.(abs.(getindex.(AF, 1)) .+ 1e-12), length(U), length(V))
af_vals = [max(af_i, -50) for af_i in af_vals]
ax3 = Axis(fig[1,3], xlabel = "u", ylabel = "v", title = "AF", aspect = DataAspect())
heatmap!(ax3, U, V, af_vals, colormap = :jet, colorrange = (-50, 0))
fig
end
