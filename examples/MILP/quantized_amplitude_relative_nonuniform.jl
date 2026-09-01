using ArraySynthesis
using ArraySynthesis: °, dB
using HiGHS

N = 24
array = symmetric_linear_array(N; d = 0.5)

mainlobe = -asin(0.15) .. asin(0.15)
sll_region = reduce(join_regions, region.(outside(mainlobe), 1°))
p = pattern(beam(0°), sidelobes(sll_region, -18dB))
obj = MinSLL(sll_region; upper_bound = -18dB)

bounds = ArraySynthesis.amplitude_upper_bounds(array, p)
@show bounds.global_bound bounds.condition

levels = [-1.0, -0.6, -0.3, 0.0, 0.3, 0.7, 1.0]
qa = ArraySynthesis.QuantizedAmplitude(levels; β = 0°, relative = true)
result = synthesize(array, p, obj, qa, MILP(big_m = bounds.global_bound), HiGHS.Optimizer; time_limit = 90.0)
@show result.status result.objective_value

result_cont = synthesize(array, p, obj, ConjugateSymmetricWeights(), SOCP(), Mosek.Optimizer)

theta_vals = -π/2:0.001:π/2
af_quantized = 20 .* log10.(max.([abs(array_factor(array, qa, result.weights, [ThetaDirection(θ)])[1]) for θ in theta_vals], 1e-12))
af_continuous = 20 .* log10.(max.([abs(array_factor(array, ConjugateSymmetricWeights(), result_cont.weights, [ThetaDirection(θ)])[1]) for θ in theta_vals], 1e-12))

using GLMakie
fig = Figure()
ax = Axis(fig[1,1], xlabel = "θ (rad)", ylabel = "|AF| (dB)")
lines!(ax, theta_vals, af_quantized, linewidth = 2, label = "Quantized")
lines!(ax, theta_vals, af_continuous, linewidth = 2, label = "Continuous")
ylims!(-60, 4); axislegend(position = :lt)
ax2 = Axis(fig[1,2])
stem!(ax2, 1:length(result.weights), real.(result.weights))
stem!(ax2, 1:length(result.weights), real.(result_cont.weights))
fig

