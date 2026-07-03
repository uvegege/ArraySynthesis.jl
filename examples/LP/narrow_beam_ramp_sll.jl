using ArraySynthesis
using ArraySynthesis: °, dB
using GLMakie
using LinearAlgebra
using HiGHS


sl_region1 = region(-90°..(-24)°, 1°)
sl_region2 = region(-4°..15°, 1°)
sl_region3 = region(15°..90°, 1°)

p = pattern(beam(-13°),
    sidelobes(sl_region1, -40dB),
    sidelobes(sl_region2, theta_ramp(-4°, -20dB, 15°, -25dB)),
    sidelobes(sl_region3, -40dB))

array = uniform_linear_array(20, d = 0.5)
obj = MaxGain()
solv = HiGHS.Optimizer
coef = ComplexWeights()
result = synthesize(array, p, obj, coef, LP(), solv)

theta_vals = -π/2:0.001:π/2
af_vals = [abs(array_factor(array, coef, result.weights, [ThetaDirection(θ)])[1]) for θ in theta_vals]
af_db = 20 .* log10.(max.(af_vals, 1e-12))

fig = Figure()
ax = Axis(fig[1,1], xlabel="θ (rad)", ylabel="|AF| (dB)")
lines!(ax, theta_vals, af_db, linewidth=2)
ylims!(-60, 10)
xlims!(-90°, 90°)
fig