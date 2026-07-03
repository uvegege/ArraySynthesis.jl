using ArraySynthesis
using ArraySynthesis: °, dB
using GLMakie
using LinearAlgebra
using HiGHS



beam_region = region(12.5°..37.5°, 1°)
sll_region1 = region(-90°..6.5°, 1°)
sll_region2 = region(43.5°..90°, 1°)

p = pattern(shaped_beam(beam_region, 1.0, ripple = -0.6dB))
obj = MinSLL(join_regions(sll_region1, sll_region2))
array = symmetric_linear_array(32, d = 0.5)

result_prog = synthesize(array, p, obj, ProgressivePhaseAmplitude(deg2rad(25)), LP(), HiGHS.Optimizer)
result_cplx = synthesize(array, p, obj, ConjugateSymmetricWeights(), LP(), HiGHS.Optimizer)

theta_vals = -π/2:0.001:π/2
coef1 = ProgressivePhaseAmplitude(deg2rad(25))
coef2 = ConjugateSymmetricWeights()
af_db1 = 20 .* log10.(max.([abs(array_factor(array, coef1, result_prog.weights, [ThetaDirection(θ)])[1]) for θ in theta_vals], 1e-12))
af_db2 = 20 .* log10.(max.([abs(array_factor(array, coef2, result_cplx.weights, [ThetaDirection(θ)])[1]) for θ in theta_vals], 1e-12))

fig = Figure()
ax = Axis(fig[1,1], xlabel="θ (rad)", ylabel="|AF| (dB)")
lines!(ax, theta_vals, af_db1, linewidth=2, label="Progressive phase")
lines!(ax, theta_vals, af_db2, linewidth=2, label="Complex weights")
ylims!(-60, 4); axislegend(position = :lt); fig
