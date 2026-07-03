using ArraySynthesis
using ArraySynthesis: °, dB
using GLMakie
using LinearAlgebra
using HiGHS
using Mosek, MosekTools


function triangular_array_example(nx::Integer, ny::Integer; dx = 0.5, dy = 0.5)
    xs = collect(-(nx ÷ 2):(nx ÷ 2)) .* dx
    ys = collect(0:ny-1) .* dy
    ys .-= maximum(ys) / 2

    p = zeros(Float64, 3, nx * ny)
    k = 1
    for iy in 1:ny
        shift = isodd(iy) ? dx / 2 : 0.0
        for ix in 1:nx
            p[1, k] = xs[ix] + shift
            p[2, k] = ys[iy]
            k += 1
        end
    end
    return ArrayGeometry(p, 2)
end

function triangle_f(center; base, height, bp = 0.0)
    u0, v0 = center
    return Polygon(
        [u0 - base / 2 - bp, u0 + base / 2 + bp, u0, u0 - base / 2 - bp],
        [v0 - height / 2 - bp / 2, v0 - height / 2 - bp / 2, v0 + height / 2 + bp, v0 - height / 2 - bp / 2],
    )
end

array = triangular_array_example(17, 17, dx = 0.5, dy = 0.5)

beam_shape = triangle_f((0.0, 0.0); base = 27°, height = 27°)
guard_shape = triangle_f((0.0, 0.0); base = 27°, height = 27°, bp = 0.2)
beam_region = region(beam_shape, step = 2°)
sl_region  = visible_region(guard_shape; step = 5°, bandpass = 0.0, filtered = false)

p = pattern(shaped_beam(beam_region, 1.0, ripple = -1.5dB),
    sidelobes(sl_region, -27dB))

obj = IterativeReweightedL1(max_iter = 15)
coef = ComplexWeights()
@time result = synthesize(array, p, obj, coef, SOCP(), Mosek.Optimizer)

active = abs.(result.weights) .> 1e-3


begin
fig = Figure()
ax = Axis(fig[1,1], xlabel="x/λ", ylabel="y/λ", title="Active: $(sum(active)) / $(length(active))")
scatter!(ax, array.positions[1, .!active], array.positions[2, .!active], color=:red, marker=:x, markersize=6)
scatter!(ax, array.positions[1,   active], array.positions[2,   active], color=:green, markersize=8)
ax2 = Axis(fig[1,2], xlabel = "u", ylabel = "v", title = "Regions", aspect = DataAspect())
scatter!(ax2, map(uv->(uv.u, uv.v), sl_region.points))
scatter!(ax2, map(uv->(uv.u, uv.v), beam_region.points))
ax3 = Axis(fig[1,3], xlabel = "u", ylabel = "v", title = "AF", aspect = DataAspect())
U = collect(-1.0:0.01:1.0)
V = collect(-1.0:0.01:1.0)
dirs = [UVDirection(u,v) for u in U, v in V]
AF = array_factor(array, coef, result.weights, dirs)
af_vals = reshape(20 .* log10.(abs.(getindex.(AF, 1)) .+ 1e-12), length(U), length(V))
af_vals = [max(af_i, -50) for af_i in af_vals]
heatmap!(ax3, U, V, af_vals, colormap = parula_cm, colorrange = (-50, 0))
ax3 = Axis3(fig[2,2])
surface!(ax3, U, V, af_vals, colormap = parula_cm)
fig
end
