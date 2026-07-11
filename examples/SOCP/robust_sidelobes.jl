using ArraySynthesis
using ArraySynthesis: °, dB
using GLMakie
using LinearAlgebra
using Random
using Clarabel

array = uniform_linear_array(16, d = 0.5)
coef = ComplexWeights()

main_lobe = ClosedInterval(-12°, 12°)
sll_region = join_regions(region.(outside(main_lobe), 1°)...)

p = pattern(
    beam(0°),
    sidelobes(sll_region, -25dB),
)

tol = robust(
    pointing_accuracy = 0.1°,
    phase_tolerance = 0.5°,
    amplitude_tolerance = 0.05dB,
    position_tolerance = (xy = 0.0005, z = 0.0),
)

result_nominal = synthesize(array, p, MinPower(), coef, SOCP(), Clarabel.Optimizer)
result_robust = synthesize(array, p, MinPower(), coef, SOCP(), Clarabel.Optimizer; robustness = tol)

theta_vals = -π/2:0.001:π/2
dirs = [ThetaDirection(θ) for θ in theta_vals]

af_nominal = abs.(array_factor(array, coef, result_nominal.weights, dirs))
af_robust = abs.(array_factor(array, coef, result_robust.weights, dirs))
af_nominal_db = 20 .* log10.(max.(af_nominal, 1e-12))
af_robust_db = 20 .* log10.(max.(af_robust, 1e-12))

certified = [
    abs(array_factor(array, coef, result_robust.weights, [dir])[1]) +
    ArraySynthesis.delta_at(tol, array, dir) * norm(result_robust.weights)
    for dir in dirs
]
certified_db = 20 .* log10.(max.(certified, 1e-12))

rng = MersenneTwister(42)
mc_trials = 120

function perturbed_pattern_db(rng, array, coef, weights, dirs, tol)
    N = size(array.positions, 2)

    Δθ = (2rand(rng) - 1) * tol.pointing[1]
    perturbed_dirs = [ThetaDirection(clamp(dir.θ + Δθ, -π / 2, π / 2)) for dir in dirs]

    amp = 1 .+ tol.amplitude .* (2rand(rng, N) .- 1)
    phase = tol.phase .* (2rand(rng, N) .- 1)
    perturbed_weights = weights .* amp .* cis.(phase)

    Δp = zeros(size(array.positions))
    for axis in 1:3
        ρ = tol.position[axis]
        ρ == 0 && continue
        Δp[axis, :] .= ρ .* (2rand(rng, N) .- 1)
    end
    perturbed_array = ArrayGeometry(array.positions .+ Δp, array.dim)

    af = abs.(array_factor(perturbed_array, coef, perturbed_weights, perturbed_dirs))
    return 20 .* log10.(max.(af, 1e-12))
end

mc_patterns_db = [
    perturbed_pattern_db(rng, array, coef, result_robust.weights, dirs, tol)
    for _ in 1:mc_trials
]

fig = Figure()
ax = Axis(fig[1, 1], xlabel = "θ (deg)", ylabel = "|AF| (dB)")
for (i, af_mc_db) in enumerate(mc_patterns_db)
    label = i == 1 ? "Monte Carlo" : nothing
    lines!(ax, theta_vals ./ °, af_mc_db, color = (:gray35, 0.16), linewidth = 1, label = label)
end
lines!(ax, theta_vals ./ °, af_nominal_db, linewidth = 2, label = "Nominal")
lines!(ax, theta_vals ./ °, af_robust_db, linewidth = 2, label = "Robust")
lines!(ax, theta_vals ./ °, certified_db, linewidth = 2, linestyle = :dash, label = "Robust bound")
hlines!(ax, [-25], linestyle = :dot, color = :gray40)
ylims!(ax, -60, 3)
xlims!(ax, -90, 90)
axislegend(ax, position = :lb)
fig
