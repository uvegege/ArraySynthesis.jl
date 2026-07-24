using ArraySynthesis
using ArraySynthesis: °, dB
using LinearAlgebra
using HiGHS

array = uniform_linear_array(22)
coef = ComplexWeights()

beam1 = -31°
beam2 = 18°
guard1 = -38°..(-24°)
guard2 = 11°..25°

sll_regions = region.(outside([guard1, guard2]), 1°)
sll_region = reduce(join_regions, sll_regions)

p = pattern(
    beam(beam1),
    beam(beam2),
)

obj = MinSLL(sll_region)
result = synthesize(array, p, obj, coef, LP(), Clarabel.Optimizer)

theta_vals = -π/2:0.001:π/2
af_vals = [abs(array_factor(array, coef, result.weights, [ThetaDirection(θ)])[1]) for θ in theta_vals]
af_db = 20 .* log10.(max.(af_vals, 1e-12))

using GLMakie
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "θ (deg)", ylabel = "|AF| (dB)")
lines!(ax, theta_vals ./ °, af_db, linewidth = 2)
vlines!(ax, [beam1, beam2] ./ °, linestyle = :dash, color = :gray40)
ylims!(-60, 5)
xlims!(-90, 90)
fig
