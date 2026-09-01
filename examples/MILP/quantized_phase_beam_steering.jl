using ArraySynthesis
using ArraySynthesis: °, dB
using HiGHS

N = 32
array = symmetric_linear_array(N; d = 0.5)

theta0 = 15°
half = asin(0.15)
mainlobe = (theta0 - half) .. (theta0 + half)
sll_region = reduce(join_regions, region.(outside(mainlobe), 1°))
p = pattern(beam(theta0), sidelobes(sll_region, -18dB))
obj = MinSLL(sll_region; upper_bound = -18dB)

M_phase = 16
phases = collect(range(0.0, 2π; length = M_phase + 1))[1:end-1]

bounds = ArraySynthesis.amplitude_upper_bounds(array, p)
qp = ArraySynthesis.QuantizedPhase(phases; free_amplitudes = true)
result = synthesize(array, p, obj, qp, MILP(big_m = bounds.global_bound), HiGHS.Optimizer; time_limit = 90.0)

result_cont = synthesize(array, p, obj, ConjugateSymmetricWeights(), LP(), HiGHS.Optimizer)

theta_vals = -π/2:0.001:π/2
af_db = 20 .* log10.(max.([abs(array_factor(array, qp, result.weights, [ThetaDirection(θ)])[1]) for θ in theta_vals], 1e-12))
af_continuous = 20 .* log10.(max.([abs(array_factor(array, ConjugateSymmetricWeights(), result_cont.weights, [ThetaDirection(θ)])[1]) for θ in theta_vals], 1e-12))

using GLMakie
fig = Figure()
ax = Axis(fig[1,1], xlabel = "θ (rad)", ylabel = "|AF| (dB)")
lines!(ax, theta_vals, af_db, linewidth = 2, label = "QuantizedPhase")
lines!(ax, theta_vals, af_continuous, linewidth = 2, label = "Continuous")

vlines!(ax, [theta0]; color = :gray, linestyle = :dash)
ylims!(-60, 4); axislegend(position = :lt)
ax2 = Axis(fig[1,2], title = "|w[n]|")
stem!(ax2, 1:length(result.weights), abs.(result.weights))
stem!(ax2, 1:length(result.weights), abs.(result_cont.weights))
ax3 = Axis(fig[1,3], title = "angle(w[n])")
stem!(ax3, 1:length(result.weights), rad2deg.(angle.(result.weights)))
stem!(ax3, 1:length(result.weights), rad2deg.(angle.(result_cont.weights)))
fig

