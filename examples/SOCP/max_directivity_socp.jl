using ArraySynthesis
using ArraySynthesis: °, dB
using GLMakie
using LinearAlgebra
using Clarabel

# Directivity maximization via SOCP: fix Re(AF(θ0)) = 1 and minimize the
# integrated power over the region, equivalent to maximizing directivity.

array = uniform_linear_array(32, d = 0.5)

main_lobe = ClosedInterval(-10°, 10°)
sll_region = join_regions(region.(outside(main_lobe), 1°)...)
integration_region = region(ClosedInterval(-90°, 90°), 1°)

p = pattern(sidelobes(sll_region, -20dB))
obj = MaxDirectivity(0°, integration_region)
coef = ComplexWeights()

result = synthesize(array, p, obj, coef, SOCP(), Clarabel.Optimizer)

theta_vals = -π/2:0.001:π/2
dirs = [ThetaDirection(θ) for θ in theta_vals]
af = abs.(array_factor(array, coef, result.weights, dirs))
#af ./= maximum(af)
af_db = 20 .* log10.(max.(af, 1e-12))

fig = Figure()
ax = Axis(fig[1, 1], xlabel = "θ (deg)", ylabel = "Normalized |AF| (dB)")
lines!(ax, theta_vals ./ °, af_db, linewidth = 2, label = "MaxDirectivity")
hlines!(ax, [-20], linestyle = :dot, color = :gray40)
ylims!(ax, -70, 2)
xlims!(ax, -90, 90)
axislegend(ax, position = :lb)
fig
