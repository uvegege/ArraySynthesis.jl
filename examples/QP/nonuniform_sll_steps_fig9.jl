using ArraySynthesis
using ArraySynthesis: °, dB
using GLMakie
using LinearAlgebra
using HiGHS

array = uniform_linear_array(21, d = 0.5)

sll_far_left = region(-90°..(-55.5°), 0.5°)
sll_near_left = region(-55°..(-10°), 0.5°)
sll_near_right = region(10°..55°, 0.5°)
sll_far_right = region(55.5°..90°, 0.5°)

p = pattern(
    sidelobes(sll_far_left, -30dB),
    sidelobes(sll_near_left, -40dB),
    sidelobes(sll_near_right, -40dB),
    sidelobes(sll_far_right, -30dB),
)

coef = ComplexWeights()
obj = IterativeFloorSynthesis(ThetaDirection(0°); max_iter = 8, tol = 5e-3)
result = synthesize(array, p, obj, coef, QP(), HiGHS.Optimizer)

theta_vals = -π/2:0.001:π/2
af_vals = [abs(array_factor(array, coef, result.weights, [ThetaDirection(θ)])[1]) for θ in theta_vals]
af_db = 20 .* log10.(max.(af_vals, 1e-12))

fig = Figure()
ax = Axis(fig[1,1], xlabel="θ (deg)", ylabel="|AF| (dB)", title="Two-level SLL")
lines!(ax, theta_vals ./ °, af_db, linewidth=2)
lines!(ax, [-90, -55], [-30, -30], linestyle=:dash, color=:black)
lines!(ax, [-55, -10], [-40, -40], linestyle=:dash, color=:black)
lines!(ax, [10, 55], [-40, -40], linestyle=:dash, color=:black)
lines!(ax, [55, 90], [-30, -30], linestyle=:dash, color=:black)
ylims!(-80, 1); xlims!(-90, 90); fig
