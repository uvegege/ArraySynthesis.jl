using ArraySynthesis
using ArraySynthesis: °, dB
using GLMakie
using LinearAlgebra
using HiGHS
using Mosek, MosekTools


array = symmetric_planar_array(14, 14, dx = 0.5, dy = 0.5)

beam_region = region(Circle(0.2, (0.0, 0.0)), step = 2°)
null_region = region(Circle(7°, (sin(96°)*cos(45°), sin(96°)*sin(45°))), step = 2°)
sl_region = visible_region(Circle(0.2, (0.0, 0.0)),
                              Circle(7°, (sin(96°)*cos(45°), sin(96°)*sin(45°))); step = 4°, bandpass = 0.1)

p = pattern(shaped_beam(beam_region, 1.0, ripple = 3dB),
    sidelobes(null_region, -45dB),
    sidelobes(sl_region, -20dB))

obj = IterativeReweightedL1(max_iter = 3)
coef = ConjugateSymmetricWeights()
result = synthesize(array, p, obj, coef, LP(), Mosek.Optimizer)

active = abs.(result.weights) .> 1e-3

begin
fig = Figure()
ax = Axis(fig[1,1], xlabel="x/λ", ylabel="y/λ", title="Active: $(sum(active)) / $(length(active))")
scatter!(ax, array.positions[1, .!active], array.positions[2, .!active], color=:red, marker=:x, markersize=6)
scatter!(ax, array.positions[1,   active], array.positions[2,   active], color=:green, markersize=8)
fig

ax2 = Axis(fig[1,2])
scatter!(ax2, map(uv->(uv.u, uv.v), sl_region.points))
scatter!(ax2, map(uv->(uv.u, uv.v), beam_region.points))
scatter!(ax2, map(uv->(uv.u, uv.v), null_region.points))

ax3 = Axis(fig[1,3])
U = collect(-1.0:0.01:1.0)
V = collect(-1.0:0.01:1.0)
dirs = [UVDirection(u,v) for u in U, v in V]
AF = array_factor(array, coef, result.weights, dirs)
af_vals = reshape(20 .* log10.(abs.(getindex.(AF, 1)) .+ 1e-12), length(U), length(V))
af_vals = [max(af_i, -45) for af_i in af_vals]
image!(ax3, af_vals, colormap = parula_cm)
ax3 = Axis3(fig[2,2])
surface!(ax3, U, V, af_vals, colormap = parula_cm)
fig

fig
end