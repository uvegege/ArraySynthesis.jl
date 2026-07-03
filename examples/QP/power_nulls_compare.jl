using ArraySynthesis
using ArraySynthesis: °, dB
using GLMakie
using LinearAlgebra
using HiGHS, Clarabel

array = uniform_linear_array(32, d = 0.5)

sll_region1 = region(-90°..(-69°), 1°)
sll_region2 = region(-71°..(-41°), 1°)
sll_region3 = region(-39°..(-30°), 1°)
sll_region4 = region(-30°..(-21°), 1°)
sll_region5 = region(-19°..(-10°), 1°)
sll_region6 = region(-10°..(-5°), 1°)
sll_region7 = region(5°..90°, 1°)

p = pattern(beam(0°),
    sidelobes(sll_region1, -20dB), sidelobes(sll_region2, -20dB),
    sidelobes(sll_region3, -20dB), sidelobes(sll_region4, -40dB),
    sidelobes(sll_region5, -40dB), sidelobes(sll_region6, -20dB),
    sidelobes(sll_region7, -20dB),
    null(-70°, level = -55dB), null(-40°, level = -55dB), null(-20°, level = -55dB))

coef = ComplexWeights()
result_qp = synthesize(array, p, MinPower(), coef, QP(), HiGHS.Optimizer)
result_lp = synthesize(array, p, Feasible(), coef, LP(), HiGHS.Optimizer)
result_socp = synthesize(array, p, MinPower(), coef, SOCP(), Clarabel.Optimizer)

theta_vals = -π/2:0.001:π/2
af_qp  = 20 .* log10.(max.([abs(array_factor(array, coef, result_qp.weights,   [ThetaDirection(θ)])[1]) for θ in theta_vals], 1e-12))
af_lp  = 20 .* log10.(max.([abs(array_factor(array, coef, result_lp.weights,   [ThetaDirection(θ)])[1]) for θ in theta_vals], 1e-12))
af_socp = 20 .* log10.(max.([abs(array_factor(array, coef, result_socp.weights, [ThetaDirection(θ)])[1]) for θ in theta_vals], 1e-12))

fig = Figure()
ax = Axis(fig[1,1], xlabel="θ (rad)", ylabel="|AF| (dB)")
lines!(ax, theta_vals, af_qp,   linewidth=2, label="QP")
lines!(ax, theta_vals, af_lp,   linewidth=2, label="LP", linestyle = :dash)
lines!(ax, theta_vals, af_socp, linewidth=2, label="SOCP", linestyle = :dot)
ylims!(-60, 1); axislegend(); fig
