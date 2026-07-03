using ArraySynthesis
using ArraySynthesis: °, dB
using GLMakie
using LinearAlgebra
using HiGHS


# Tseng-Griffiths Fig. 3-style: 20-element half-wavelength ULA, look
# direction 50°, and equiripple sidelobe target -40 dB. This uses the same
# angular specification as focused_nonuniform.jl, but with uniform spacing.
array = uniform_linear_array(20, d = 0.5)

sll_region1 = region(-90°..34°, 0.5°)
sll_region2 = region(74°..90°, 0.5°)
p = pattern(sidelobes(sll_region1, -40dB), sidelobes(sll_region2, -40dB))

coef = ComplexWeights()
obj = IterativeFloorSynthesis(ThetaDirection(50°); max_iter = 8)
result = synthesize(array, p, obj, coef, QP(), HiGHS.Optimizer)

theta_vals = -π/2:0.001:π/2
af_vals = [abs(array_factor(array, coef, result.weights, [ThetaDirection(θ)])[1]) for θ in theta_vals]
af_db = 20 .* log10.(max.(af_vals, 1e-12))

fig = Figure()
ax = Axis(fig[1,1], xlabel="θ (deg)", ylabel="|AF| (dB)", title="20-element ULA, look 50°")
lines!(ax, theta_vals ./ °, af_db, linewidth=2)
hlines!(ax, [-40], linestyle=:dash, color=:black)
ylims!(-80, 1); xlims!(-90, 90); fig
