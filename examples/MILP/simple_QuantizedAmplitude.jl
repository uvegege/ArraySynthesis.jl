cd(raw"C:\MisProyecto\ConvexAntennaArray\Julia\ArraySynthesis-quantized")
using Pkg
Pkg.activate(".")
using ArraySynthesis
using ArraySynthesis: °, dB
Pkg.activate(; temp = true)
Pkg.add(["HiGHS", "Clarabel", "MosekTools", "Mosek"])
using HiGHS, Clarabel, MosekTools, Mosek
Pkg.activate(".")

N = 20
array = symmetric_linear_array(N, d = 0.5)

mainlobe = -asin(0.22) .. asin(0.22)
sll_region = reduce(join_regions, region.(outside(mainlobe), 1°))

p = pattern(beam(0°), sidelobes(sll_region, -18dB))
obj = MinSLL(sll_region; upper_bound = -18dB)

#levels = [0.03, 0.05, 0.07]
#qw = ArraySynthesis.QuantizedAmplitude(levels; β = 0°, relative = false)
N_values = 4
levels =  range(0.0, 1.0, N_values)
qw = ArraySynthesis.QuantizedAmplitude(levels; β = 0°, relative = true)
milp_result = synthesize(array, p, obj, qw, MILP(), HiGHS.Optimizer)
theta_vals = -π/2:0.0005:π/2
af_vals = [abs(array_factor(array, qw, milp_result.weights, [ThetaDirection(θ)])[1]) for θ in theta_vals]
af_db = 20 .* log10.(max.(af_vals, 1e-12))

using GLMakie
fig = Figure()
ax = Axis(fig[1,1], xlabel="sin(θ)", ylabel="|AF| (dB)")
lines!(ax, sin.(theta_vals), af_db, linewidth=2)
xlims!(-1, 1)
ylims!(-80, 5)
ax2 = Axis(fig[1,2])
stem!(ax2, 1:length(milp_result.weights), real.(milp_result.weights))
fig



# Flat top case
beam_region = region(12.5°..37.5°, 1°)
sll_region1 = region(-90°..6.5°, 1°)
sll_region2 = region(43.5°..90°, 1°)

p = pattern(shaped_beam(beam_region, 1.0, ripple = -1.0dB))
obj = MinSLL(join_regions(sll_region1, sll_region2))
array = symmetric_linear_array(32, d = 0.5)

#=
N_values = 32
levels = collect(range(-1.0, 1.0; length = N_values+1))
qa = ArraySynthesis.QuantizedAmplitude(levels, β = deg2rad(25), relative = true)
result = synthesize(array, p, obj, qa, MILP(), Mosek.Optimizer, time_limit = 60)

using GLMakie
theta_vals = -π/2:0.001:π/2
af_db1 = 20 .* log10.(max.([abs(array_factor(array, qa, result.weights, [ThetaDirection(θ)])[1]) for θ in theta_vals], 1e-12))

fig = Figure()
ax = Axis(fig[1,1], xlabel="θ (rad)", ylabel="|AF| (dB)")
lines!(ax, theta_vals, af_db1, linewidth=2, label="Quantized Amplitudes")
ylims!(-60, 4); axislegend(position = :lt)
ax2 = Axis(fig[1,2])
stem!(ax2, 1:length(result.weights), real.(result.weights))
fig
=#

# `MinSLL` can be slow compared to `Feasible()`
N_values = 32
levels = collect(range(-1.0, 1.0; length = N_values+1))
qa = ArraySynthesis.QuantizedAmplitude(levels, β = deg2rad(25), relative = true)
p2 = pattern(shaped_beam(beam_region, 1.0, ripple = -1.0dB), 
             sidelobes(join_regions(sll_region1, sll_region2), -20dB))

result2 = synthesize(array, p2, Feasible(), qa, MILP(big_m = 100), Mosek.Optimizer, time_limit = 60)
@show ArraySynthesis.JuMP.solve_time(result2.model)

theta_vals = -π/2:0.001:π/2
af_db1 = 20 .* log10.(max.([abs(array_factor(array, qa, result2.weights, [ThetaDirection(θ)])[1]) for θ in theta_vals], 1e-12))

using GLMakie
fig = Figure()
ax = Axis(fig[1,1], xlabel="θ (rad)", ylabel="|AF| (dB)")
lines!(ax, theta_vals, af_db1, linewidth=2, label="Quantized Amplitudes")
ylims!(-60, 4); axislegend(position = :lt)
ax2 = Axis(fig[1,2])
stem!(ax2, 1:length(result2.weights), real.(result2.weights))
fig
