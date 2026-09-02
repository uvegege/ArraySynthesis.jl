using ArraySynthesis
using ArraySynthesis: °, dB
using MosekTools, Mosek

N = 32
array = materialize(symmetric_linear_array(N; d = 0.5))

beam_region = region(12.5° .. 37.5°, 1°)
sll_region1 = region(-90° .. 6.5°, 1°)
sll_region2 = region(43.5° .. 90°, 1°)
sll_region = join_regions(sll_region1, sll_region2)
p = pattern(shaped_beam(beam_region, 1.0; ripple = -1.0dB), sidelobes(sll_region, -18dB))
obj = MinSLL(sll_region; upper_bound = -18dB)

result_cont = synthesize(array, p, obj, ComplexWeights(), SOCP(), Mosek.Optimizer)
@show extrema(real.(result_cont.weights)) extrema(imag.(result_cont.weights))

L = 32
re_levels = collect(range(-0.25, 0.25; length = L))
im_levels = collect(range(-0.25, 0.25; length = L))
C = ComplexF64[complex(re, im) for re in re_levels for im in im_levels]
qw = ArraySynthesis.QuantizedWeights(C; free_scale = false)

result = synthesize(array, p, obj, qw, MILP(), Mosek.Optimizer; time_limit = 120.0)

theta_vals = -π/2:0.001:π/2
af_quantized = 20 .* log10.(max.([abs(array_factor(array, qw, result.weights, [ThetaDirection(θ)])[1]) for θ in theta_vals], 1e-12))
af_continuous = 20 .* log10.(max.([abs(array_factor(array, ComplexWeights(), result_cont.weights, [ThetaDirection(θ)])[1]) for θ in theta_vals], 1e-12))

using GLMakie
fig = Figure()
ax = Axis(fig[1,1], xlabel = "θ (rad)", ylabel = "|AF| (dB)")
lines!(ax, theta_vals, af_quantized, linewidth = 2, label = "QuantizedWeights")
lines!(ax, theta_vals, af_continuous, linewidth = 2, label = "Continuous")
ylims!(-60, 4); axislegend(position = :lt)
ax2 = Axis(fig[1,2], title = "|w[n]|")
stem!(ax2, 1:length(result.weights), abs.(result.weights))
stem!(ax2, 1:length(result_cont.weights), abs.(result_cont.weights))
ax3 = Axis(fig[1,3], title = "angle(w[n])")
stem!(ax3, 1:length(result.weights), rad2deg.(angle.(result.weights)))
stem!(ax3, 1:length(result_cont.weights), rad2deg.(angle.(result_cont.weights)))
fig

