# # Narrow beam with a ramp sidelobe mask
#
# This example uses `MaxGain` with explicit sidelobe constraints. The objective
# maximizes the response in the look direction, while the mask controls what is
# allowed elsewhere.

using ArraySynthesis
using ArraySynthesis: °, dB
using GLMakie
using HiGHS

array = uniform_linear_array(20, d = 0.5)
coef = ComplexWeights()

# The main beam is a single look direction. The middle sidelobe region uses a
# functional mask, so the allowed level changes with angle.
sl_region1 = region(ClosedInterval(-90°, -24°), 1°)
sl_region2 = region(ClosedInterval(-4°, 15°), 1°)
sl_region3 = region(ClosedInterval(15°, 90°), 1°)

p = pattern(
    beam(-13°),
    sidelobes(sl_region1, -40dB),
    sidelobes(sl_region2, theta_ramp(-4°, -20dB, 15°, -25dB)),
    sidelobes(sl_region3, -40dB),
)

# With no explicit direction, `MaxGain()` uses the beam direction from `p`.
result = synthesize(array, p, MaxGain(), coef, LP(), HiGHS.Optimizer)

theta_vals = -π/2:0.001:π/2
af_vals = [abs(array_factor(array, coef, result.weights, [ThetaDirection(θ)])[1]) for θ in theta_vals]
af_db = 20 .* log10.(max.(af_vals, 1e-12))

fig = Figure()
ax = Axis(fig[1, 1], xlabel = "θ (deg)", ylabel = "|AF| (dB)")
lines!(ax, theta_vals ./ °, af_db, linewidth = 2)
ylims!(ax, -60, 10)
xlims!(ax, -90, 90)
fig

# The documentation includes a precomputed result image, so this example is
# shown without being executed by Documenter.
#
# ![Narrow beam with ramp sidelobe mask](../assets/narrow_beam_ramp_sll.png)
