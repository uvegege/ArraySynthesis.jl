using ArraySynthesis
using ArraySynthesis: °, dB
using GLMakie
using LinearAlgebra
using HiGHS

posx = [0, 0.4618, 0.8425, 1.2126, 1.6328, 1.9412, 2.5037, 3.0292, 3.4567, 3.9539,
        4.2221, 4.5851, 5.2431, 5.6073, 6.2850, 6.5660, 7.1697, 7.5401, 8.1355, 8.4386, 8.8206]
positions = zeros(3, length(posx)); positions[1, :] .= posx
array = ArrayGeometry(positions, 1)

sll_region1 = region(-90°..(-10°), 0.5°)
sll_region2 = region(10°..90°, 0.5°)
p = pattern(sidelobes(sll_region1, -40dB), sidelobes(sll_region2, -40dB))

coef = ComplexWeights()
obj = IterativeFloorSynthesis(ThetaDirection(0°); max_iter = 8)
result = synthesize(array, p, obj, coef, QP(), HiGHS.Optimizer)

theta_vals = -π/2:0.001:π/2
af_vals = [abs(array_factor(array, coef, result.weights, [ThetaDirection(θ)])[1]) for θ in theta_vals]
af_db = 20 .* log10.(max.(af_vals, 1e-12))

fig = Figure()
ax = Axis(fig[1,1], xlabel="θ (deg)", ylabel="|AF| (dB)", title="Random array, look 0°")
lines!(ax, theta_vals ./ °, af_db, linewidth=2)
hlines!(ax, [-40], linestyle=:dash, color=:black)
ylims!(-80, 1); xlims!(-90, 90); fig
