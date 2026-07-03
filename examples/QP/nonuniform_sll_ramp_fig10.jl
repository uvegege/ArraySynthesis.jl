using ArraySynthesis
using ArraySynthesis: °, dB
using GLMakie
using LinearAlgebra
using HiGHS

array = uniform_linear_array(21, d = 0.5)

ramp_level(p) = begin
    θdeg = p.θ / °
    θdeg < 40 ? ((-40 - (-30)) / (40 - (-90))) * (θdeg - (-90)) - 30 : -30
end

sll_region1 = region(-90°..40°, 0.5°)
sll_region2 = region(60°..90°, 0.5°)
p = pattern(
    sidelobes(sll_region1, p -> ramp_level(p) * dB),
    sidelobes(sll_region2, -30dB),
)

coef = ComplexWeights()
obj = IterativeFloorSynthesis(ThetaDirection(50°); max_iter = 8)
result = synthesize(array, p, obj, coef, QP(), HiGHS.Optimizer)

theta_vals = -π/2:0.001:π/2
af_vals = [abs(array_factor(array, coef, result.weights, [ThetaDirection(θ)])[1]) for θ in theta_vals]
af_db = 20 .* log10.(max.(af_vals, 1e-12))

fig = Figure()
ax = Axis(fig[1,1], xlabel="θ (deg)", ylabel="|AF| (dB)", title="Ramp SLL, look 50°")
lines!(ax, theta_vals ./ °, af_db, linewidth=2)
spec_left = collect(-90:1:40)
lines!(ax, spec_left, [ramp_level(ThetaDirection(θ * °)) for θ in spec_left], linestyle=:dash, color=:black)
lines!(ax, [60, 90], [-30, -30], linestyle=:dash, color=:black)
ylims!(-80, 1); xlims!(-90, 90); fig
