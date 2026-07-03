using ArraySynthesis
using ArraySynthesis: °, dB
using GLMakie
using LinearAlgebra
using Mosek
using MosekTools
using HiGHS

array = planar_array(20, 20, dx = 0.5, dy = 0.5)

beam1_dir = uv(37°, 37°)
beam2_dir = uv(-37°, 37°)
beam_shape1 = Rectangle(9°, 9°, (beam1_dir.u, beam1_dir.v))
beam_shape2 = Rectangle(9°, 9°, (beam2_dir.u, beam2_dir.v))
beam_region1 = region(beam_shape1, step = 2°)
beam_region2 = region(beam_shape2, step = 2°)
sl_region = visible_region(beam_shape1, beam_shape2; step = 4°, bandpass = 0.1)

p = pattern(shaped_beam(beam_region1, 1.0, ripple = -1.5dB),
    shaped_beam(beam_region2, 1.0, ripple = -1.5dB))

obj = MinSLL([sl_region])
coef = ComplexWeights()
@time result = synthesize(array, p, obj, coef, LP(), Mosek.Optimizer)
#@time result = synthesize(array, p, obj, coef, LP(), HiGHS.Optimizer; solver_options = Dict("solver" => "ipm"))

U = V = collect(-1.0:0.02:1.0)
dirs = [UVDirection(u, v) for u in U, v in V]
AF = array_factor(array, coef, result.weights, dirs)
af_vals = reshape(20 .* log10.(abs.(getindex.(AF, 1)) .+ 1e-12), length(U), length(V))

fig = Figure()
ax = Axis(fig[1,1], xlabel="u", ylabel="v")
image!(ax, (-1,1), (-1,1), af_vals, colorrange=(-50, 0), colormap = parula_cm)
Colorbar(fig[1,3], colorrange=(-60, 0)); fig
ax2 = Axis3(fig[1,2])
surface!(ax2, U, V, af_vals, colormap = parula_cm)
fig
