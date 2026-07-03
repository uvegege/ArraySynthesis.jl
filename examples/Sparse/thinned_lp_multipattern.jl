using ArraySynthesis
using ArraySynthesis: °, dB
using GLMakie
using LinearAlgebra
using HiGHS
using Mosek, MosekTools

array = planar_array(14, 14, dx = 0.5, dy = 0.5)

beam_region1 = region(Circle(0.2, (0.0, 0.0)), step = 2°)
sl_region1 = visible_region(Circle(0.4, (0.0, 0.0)); step = 4°, bandpass = 0.0, filtered = false)
p1 = pattern(
    shaped_beam(beam_region1, 1.0, ripple = -1dB),
    sidelobes(sl_region1, -25.85dB),
)

beam_shape2 = rhombus((0.2, 0.2), 0.2)
guard_shape2 = rhombus((0.2, 0.2), 0.4)
beam_region2 = region(beam_shape2, step = 2°)
sl_region2 = visible_region(guard_shape2, Circle(0.1, (-0.5, -0.5)); step = 4°, bandpass = 0.0, filtered = false)
null_region2 = region(Circle(0.1, (-0.5, -0.5)), step = 2°)

p2 = pattern(
    shaped_beam(beam_region2, 1.0, ripple = -1dB),
    sidelobes(sl_region2, -24.30dB),
    sidelobes(null_region2, -50dB),
)

obj = MultiPatternReweightedL1(max_iter = 15)
coef = ComplexWeights()
result = synthesize(array, [p1, p2], obj, coef, LP(), Mosek.Optimizer)

# An element is active if it is used by any pattern
max_activity = [maximum(abs(result.weights[k][n]) for k in eachindex(result.weights)) for n in eachindex(result.weights[1])]
active = max_activity .> 1e-5

fig = Figure()
ax = Axis(fig[1,1], xlabel="x/λ", ylabel="y/λ", title="Active: $(sum(active)) / $(length(active))")
scatter!(ax, array.positions[1, .!active], array.positions[2, .!active], color=:red, marker=:x, markersize=6)
scatter!(ax, array.positions[1,   active], array.positions[2,   active], color=:green, markersize=8)
fig

fig_regions = Figure(size = (900, 400))
region_sets = [
    ("pattern 1 regions", beam_region1, sl_region1, nothing),
    ("pattern 2 regions", beam_region2, sl_region2, null_region2),
]

for k in 1:2
    title, beam_region, sl_region, null_region = region_sets[k]
    axr = Axis(fig_regions[1, k], xlabel = "u", ylabel = "v", title = title, aspect = DataAspect())

    sl_region !== nothing &&
        scatter!(axr, map(uv -> (uv.u, uv.v), sl_region.points), color = (:red, 0.35), markersize = 3)
    null_region !== nothing &&
        scatter!(axr, map(uv -> (uv.u, uv.v), null_region.points), color = (:orange, 0.65), markersize = 5)
    scatter!(axr, map(uv -> (uv.u, uv.v), beam_region.points), color = :green, markersize = 5)

    xlims!(axr, -1, 1)
    ylims!(axr, -1, 1)
end
fig_regions

U = V = collect(-1.0:0.02:1.0)
dirs = [UVDirection(u, v) for u in U, v in V]
inside_visible = [u^2 + v^2 <= 1.0 for u in U, v in V]

function pattern_db(array, coef, weights, dirs, nU, nV)
    af = array_factor(array, coef, weights, dirs)
    values = reshape(20 .* log10.(abs.(af) .+ 1e-12), nU, nV)
    return values .- maximum(values)
end

fig_af = Figure(size = (900, 420))
for k in 1:2
    axa = Axis(fig_af[1, k], xlabel = "u", ylabel = "v", title = "pattern $k AF", aspect = DataAspect())
    af_db = pattern_db(array, coef, result.weights[k], dirs, length(U), length(V))
    af_db[.!inside_visible] .= NaN
    heatmap!(axa, U, V, af_db, colorrange = (-40, 0), colormap = :jet)
    xlims!(axa, -1, 1)
    ylims!(axa, -1, 1)
end
Colorbar(fig_af[1, 3], colorrange = (-40, 0), label = "|AF| (dB)", colormap = :jet)
fig_af
