using ArraySynthesis
using ArraySynthesis: °, dB
using HiGHS

N = 17
array = symmetric_linear_array(N; d = 0.5)

mainlobe = -asin(0.15) .. asin(0.15)
sll_region = reduce(join_regions, region.(outside(mainlobe), 1°))
p = pattern(beam(0°))#, sidelobes(sll_region, -18dB))

levels = collect(range(-1.0, 1.0; length = 8))
qa = ArraySynthesis.QuantizedAmplitude(levels; β = 0°, relative = true)
result = synthesize(array, p, MinSLL(sll_region), qa, MILP(big_m = 1.0), HiGHS.Optimizer; time_limit = 60.0)

peak_sll_db = 20 * log10(maximum(abs.(array_factor(array, qa, result.weights, sll_region.points))))

theta_vals = -π/2:0.001:π/2
af_db = 20 .* log10.(max.([abs(array_factor(array, qa, result.weights, [ThetaDirection(θ)])[1]) for θ in theta_vals], 1e-12))

using GLMakie
fig = Figure()
ax = Axis(fig[1,1], xlabel = "θ (rad)", ylabel = "|AF| (dB)")
lines!(ax, theta_vals, af_db, linewidth = 2, label = "QuantizedAmplitude (relative, uniform)")
ylims!(-60, 4); axislegend(position = :lt)
ax2 = Axis(fig[1,2])
stem!(ax2, 1:length(result.weights), abs.(result.weights))
fig

