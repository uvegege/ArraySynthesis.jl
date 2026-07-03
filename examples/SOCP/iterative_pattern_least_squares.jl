using ArraySynthesis
using ArraySynthesis: °, dB
using GLMakie
using LinearAlgebra
using Clarabel

array = uniform_linear_array(28, d = 0.5)

beam_region = region(-25°..25°, 1°)
sll_region1 = region(-90°..(-32°), 1°)
sll_region2 = region(32°..90°, 1°)

target = [0.55 + 0.45 * cos(π * p.θ / 25°) for p in beam_region.points]

p = pattern(
    sidelobes(sll_region1, -25dB),
    sidelobes(sll_region2, -25dB),
)
obj = IterativePatternLeastSquares(beam_region, target; max_iter = 12)
coef = ComplexWeights()

result = synthesize(array, p, obj, coef, SOCP(), Clarabel.Optimizer)

theta_vals = -π/2:0.002:π/2
af_vals = [abs(array_factor(array, coef, result.weights, [ThetaDirection(θ)])[1]) for θ in theta_vals]
af_db = 20 .* log10.(max.(af_vals, 1e-12))
target_db = 20 .* log10.(max.(target, 1e-12))

fig = Figure()
ax = Axis(fig[1, 1], xlabel = "θ (deg)", ylabel = "|AF| (dB)")
lines!(ax, theta_vals ./ °, af_db, linewidth = 2, label = "Iterative LS")
lines!(ax, [p.θ / ° for p in beam_region.points], target_db, linestyle = :dash, color = :black, label = "Target")
hlines!(ax, [-25], linestyle = :dot, color = :gray40)
ylims!(-50, 2)
xlims!(-90, 90)
axislegend(position = :lb)
fig
