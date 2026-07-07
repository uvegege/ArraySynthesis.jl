using ArraySynthesis
using ArraySynthesis: °, dB
using LinearAlgebra
using HiGHS
using ArraySynthesis: ..

# Ng et al., "Flexible Array Synthesis via Quadratic Programming",
# desired-pattern approximation style example.

array = uniform_linear_array(32, d = 0.5)

fit_region = region(-90°..90°, 1°)
target = [
    abs(p.θ) <= 12° ? 1.0 :
    abs(p.θ) <= 18° ? (18° - abs(p.θ)) / 6° :
    10^(-35 / 20)
    for p in fit_region.points
]

p = pattern(beam(0°))
obj = MinFieldError(fit_region, target)
coef = ComplexWeights()

result = synthesize(array, p, obj, coef, QP(), HiGHS.Optimizer)

theta_vals = -π/2:0.001:π/2
dirs = [ThetaDirection(θ) for θ in theta_vals]
af = abs.(array_factor(array, coef, result.weights, dirs))
af ./= maximum(af)
af_db = 20 .* log10.(max.(af, 1e-12))
target_db = 20 .* log10.(max.(target, 1e-12))

using GLMakie
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "θ (deg)", ylabel = "Normalized |AF| (dB)")
lines!(ax, theta_vals ./ °, af_db, linewidth = 2, label = "MinFieldError")
lines!(ax, [p.θ / ° for p in fit_region.points], target_db, linestyle = :dash, color = :black, label = "Target")
ylims!(ax, -70, 2)
xlims!(ax, -90, 90)
axislegend(ax, position = :lb)
fig
