using ArraySynthesis
using ArraySynthesis: °, dB
using GLMakie
using LinearAlgebra
using HiGHS

# Ng et al., "Flexible Array Synthesis via Quadratic Programming",
# minimum-average sidelobe style example.

array = symmetric_linear_array(16, d = 0.5)

main_lobe = -12°..12°
sll_region = join_regions(region.(outside(main_lobe), 1°)...)

p = pattern(beam(0°))
obj = MinIntegratedPower(sll_region)
coef = ConjugateSymmetricWeights()

result = synthesize(array, p, obj, coef, QP(), HiGHS.Optimizer)

theta_vals = -π/2:0.001:π/2
dirs = [ThetaDirection(θ) for θ in theta_vals]
af = abs.(array_factor(array, coef, result.weights, dirs))
#af ./= maximum(af)
af_db = 20 .* log10.(max.(af, 1e-12))

fig = Figure()
ax = Axis(fig[1, 1], xlabel = "θ (deg)", ylabel = "Normalized |AF| (dB)")
lines!(ax, theta_vals ./ °, af_db, linewidth = 2, label = "MinIntegratedPower")
hlines!(ax, [-30], linestyle = :dot, color = :gray40)
ylims!(ax, -70, 2)
xlims!(ax, -90, 90)
axislegend(ax, position = :lb)
fig
