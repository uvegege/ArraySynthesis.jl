using ArraySynthesis
using ArraySynthesis: °, dB
using GLMakie
using LinearAlgebra
using HiGHS
using Mosek, MosekTools


# Triangular-grid planar array from:
# Synthesis of Sparse or Thinned Linear and Planar Arrays Generating
# Reconfigurable Multiple Real Patterns by Iterative Linear Programming.
#
# Article example 3: 247 potential positions, dx = lambda / sqrt(3),
# dy = lambda / 2, four reconfigurable real patterns, ripple <= 1 dB
# and SLL <= -20 dB.

const ripple = -1.2dB
const sll = -20dB

function centered_triangular_array(Nx::Integer, Ny::Integer; dx = 1 / sqrt(3), dy = 0.5)
    isodd(Nx) && isodd(Ny) || error("Centered triangular symmetry expects odd Nx and Ny.")
    p = zeros(Float64, 3, Nx * Ny)
    rows = collect(-(Ny ÷ 2):(Ny ÷ 2))
    cols = collect(-(Nx ÷ 2):(Nx ÷ 2))

    k = 1
    for r in rows
        shift = isodd(r) ? 0.5 * sign(r) : 0.0
        for c in cols
            p[1, k] = (c + shift) * dx
            p[2, k] = r * dy
            k += 1
        end
    end

    return ArrayGeometry(p, 2)
end

array_full = centered_triangular_array(13, 19)
array = SymmetricArray(symmetrize(array_full.positions)...)

focused_main = uv(0°, 0°)
focused_guard = Circle(0.17, (0.0, 0.0))
focused_sl = visible_region(focused_guard; step = 4°, bandpass = 0.0, filtered = false)
p1 = pattern(beam(focused_main), sidelobes(focused_sl, sll))

circular_main = Circle(0.2, (0.0, 0.0))
circular_guard = Circle(0.4, (0.0, 0.0))
circular_region = region(circular_main; step = 2°)
circular_sl = visible_region(circular_guard; step = 4°, bandpass = 0.0, filtered = false)
p2 = pattern(shaped_beam(circular_region, 1.0, ripple = ripple), sidelobes(circular_sl, sll))

moon_main = Moonlike(0.35, 0.5, (0.0, 0.0), (-0.35, 0.0))
moon_guard = Moonlike(0.5, 0.5, (0.0, 0.0), (-0.5, 0.0))
moon_region = region(moon_main; step = 2°)
moon_sl = visible_region(moon_guard; step = 4°, bandpass = 0.0, filtered = false)
p3 = pattern(shaped_beam(moon_region, 1.0, ripple = ripple), sidelobes(moon_sl, sll))

rect_main = Rectangle(0.25, 0.1, (0.0, 0.0))
rect_guard = Rectangle(0.4, 0.25, (0.0, 0.0))
rect_region = region(rect_main; step = 2°)
rect_sl = visible_region(rect_guard; step = 4°, bandpass = 0.0, filtered = false)
p4 = pattern(shaped_beam(rect_region, 1.0, ripple = ripple), sidelobes(rect_sl, sll))

patterns = [p1, p2, p3, p4]
obj = MultiPatternReweightedL1(max_iter = 9)
coef = ConjugateSymmetricWeights()
result = synthesize(array, patterns, obj, coef, LP(polygon_faces = 8), Mosek.Optimizer)

max_activity = [maximum(abs(result.weights[k][n]) for k in eachindex(result.weights)) for n in eachindex(result.weights[1])]
active_representatives = max_activity .> 1e-3

function expand_symmetric_activity(array::SymmetricArray, active_representatives)
    active = Bool[]
    for (i, p) in enumerate(eachcol(array.positions))
        push!(active, active_representatives[i])
        is_origin(p) || push!(active, active_representatives[i])
    end
    return active
end

active = expand_symmetric_activity(array, active_representatives)
full_array = materialize(array)

println("Potential positions: ", size(full_array.positions, 2))
println("Selected antennas: ", count(active), " (", round(100 * count(active) / length(active), digits = 1), "%)")
println("Converged: ", result.converged, " after ", result.iterations, " iterations")

U = V = collect(-1.0:0.01:1.0)
dirs = [UVDirection(u, v) for u in U, v in V]
inside_visible = [u^2 + v^2 <= 1.0 for u in U, v in V]

function pattern_db(array, coef, weights, dirs, nU, nV)
    af = array_factor(array, coef, weights, dirs)
    values = reshape(20 .* log10.(abs.(af) .+ 1e-12), nU, nV)
    return values .- maximum(values)
end

fig = Figure(size = (1200, 760))
ax0 = Axis(fig[1:2, 1], xlabel = "x/lambda", ylabel = "y/lambda",
    title = "Selected: $(count(active)) / $(length(active))")
scatter!(ax0, full_array.positions[1, .!active], full_array.positions[2, .!active],
    color = :red, marker = :x, markersize = 5)
scatter!(ax0, full_array.positions[1, active], full_array.positions[2, active],
    color = :green, markersize = 7)

titles = ["focused", "circular", "moon-like", "rectangular"]
for k in 1:4
    row = k <= 2 ? 1 : 2
    col = k <= 2 ? k + 1 : k - 1
    ax = Axis(fig[row, col], xlabel = "u", ylabel = "v", title = titles[k])
    af_db = pattern_db(array, coef, result.weights[k], dirs, length(U), length(V))
    af_db[.!inside_visible] .= NaN
    heatmap!(ax, U, V, af_db, colorrange = (-40, 0), colormap = parula_cm)
end

Colorbar(fig[1:2, 4], colorrange = (-40, 0), label = "|AF| (dB)", colormap = parula_cm)
fig


begin
  fig_regions = Figure(size = (900, 700))

  region_sets = [
      ("focused", nothing, focused_sl),
      ("circular", circular_region, circular_sl),
      ("moon-like", moon_region, moon_sl),
      ("rectangular", rect_region, rect_sl),
  ]

  for k in 1:4
      row = k <= 2 ? 1 : 2
      col = k <= 2 ? k : k - 2
      ax = Axis(fig_regions[row, col], xlabel = "u", ylabel = "v", title = region_sets[k][1])

      beam_region = region_sets[k][2]
      sl_region = region_sets[k][3]

      scatter!(ax, map(uv -> (uv.u, uv.v), sl_region.points), markersize = 3, color = (:red, 0.35))

      if beam_region !== nothing
          scatter!(ax, map(uv -> (uv.u, uv.v), beam_region.points), markersize = 5, color = :green)
      else
          scatter!(ax, [(focused_main.u, focused_main.v)], markersize = 10, color = :green)
      end

      xlims!(ax, -1, 1)
      ylims!(ax, -1, 1)
  end
  fig_regions
end
