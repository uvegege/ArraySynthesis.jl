using ArraySynthesis
using ArraySynthesis: °, dB
MosekTools, Mosek

N = 32
array = symmetric_linear_array(N; d = 0.5)

beam_region = region(12.5° .. 37.5°, 1°)
sll_region1 = region(-90° .. 6.5°, 1°)
sll_region2 = region(43.5° .. 90°, 1°)
sll_region = join_regions(sll_region1, sll_region2)
p = pattern(shaped_beam(beam_region, 1.0; ripple = -1.0dB), sidelobes(sll_region, -18dB))
obj = MinSLL(sll_region; upper_bound = -18dB)

result_cont = synthesize(array, p, obj, ConjugateSymmetricWeights(), LP(), Mosek.Optimizer)
taper = abs.(result_cont.weights)

M_phase = 16
phases = collect(range(0.0, 2π; length = M_phase + 1))[1:end-1]
qp1 = ArraySynthesis.QuantizedPhase(phases, free_amplitudes = true)
result_free = synthesize(array, p, obj, qp1, MILP(), Mosek.Optimizer; time_limit = 300.0)
qp2 = ArraySynthesis.QuantizedPhase(phases, amplitude = taper, free_amplitudes = false)
result_taper = synthesize(array, p, obj, qp2, MILP(), Mosek.Optimizer; time_limit = 300.0)


theta_vals = -π/2:0.001:π/2
af_quantized1 = 20 .* log10.(max.([abs(array_factor(array, qp, result_free.weights, [ThetaDirection(θ)])[1]) for θ in theta_vals], 1e-12))
af_quantized2 = 20 .* log10.(max.([abs(array_factor(array, qp, result_taper.weights, [ThetaDirection(θ)])[1]) for θ in theta_vals], 1e-12))
af_continuous = 20 .* log10.(max.([abs(array_factor(array, ConjugateSymmetricWeights(), result_cont.weights, [ThetaDirection(θ)])[1]) for θ in theta_vals], 1e-12))

using GLMakie
fig = Figure()
ax = Axis(fig[1,1], xlabel = "θ (rad)", ylabel = "|AF| (dB)")
lines!(ax, theta_vals, af_quantized1, linewidth = 2, label = "QuantizedPhase (free)")
lines!(ax, theta_vals, af_quantized2, linewidth = 2, label = "QuantizedPhase (taper)")
lines!(ax, theta_vals, af_continuous, linewidth = 2, label = "Continuous")
ylims!(-60, 4); axislegend(position = :lt)
ax2 = Axis(fig[1,2], title = "|w[n]|")
stem!(ax2, 1:length(result.weights), abs.(result_free.weights))
stem!(ax2, 1:length(result.weights), abs.(result_taper.weights))
stem!(ax2, 1:length(result_cont.weights), abs.(result_cont.weights))
ax3 = Axis(fig[1,3], title = "angle(w[n])")
stem!(ax3, 1:length(result.weights), rad2deg.(angle.(result_free.weights)))
stem!(ax3, 1:length(result.weights), rad2deg.(angle.(result_taper.weights)))
stem!(ax3, 1:length(result_cont.weights), rad2deg.(angle.(result_cont.weights)))
fig
