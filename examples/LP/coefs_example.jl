using ArraySynthesis
using ArraySynthesis: °, dB
using GLMakie
using LinearAlgebra
using HiGHS

sll_region1 = region(-90°..(-40°), 0.2°)
sll_region2 = region(-40°..(-15)°, 0.2°)
sll_region3 = region(15°..90°, 0.2°)

p = pattern(beam(0°), 
    sidelobes(sll_region1, -35dB), 
    sidelobes(sll_region2, -50dB), 
    sidelobes(sll_region3, -45dB))

posx = [0.3807, 0.9685, 1.3071, 1.9853, 2.2606, 2.8423, 3.2654, 3.7824, 4.5245, -0.3807, -0.9685, -1.3071, -1.9853, -2.2606, -2.8423, -3.2654, -3.7824, -4.5245]
positions = zeros(3, length(posx))
positions[1, :] .= posx
array = SymmetricArray(symmetrize(positions)...)

obj = MinSLL([sll_region1, sll_region2, sll_region3])
coefs = ProgressivePhaseAmplitude(0)
resultado1 = synthesize(array, p, obj, coefs, LP(), HiGHS.Optimizer)

af_vals1 = [abs(array_factor(array, coefs, resultado1.weights, [ThetaDirection(θ)])[1]) for θ in theta_vals]
af_db1 = 20 .* log10.(max.(af_vals1, 1e-12))

coefs2 = ConjugateSymmetricWeights()
resultado2 = synthesize(array, p, obj, coefs2, LP(), HiGHS.Optimizer)
af_vals2 = [abs(array_factor(array, coefs2, resultado2.weights, [ThetaDirection(θ)])[1]) for θ in theta_vals]
af_db2 = 20 .* log10.(max.(af_vals2, 1e-12))

fig = Figure()
ax = Axis(fig[1,1], xlabel="θ (º)", ylabel="|AF| (dB)")
lines!(ax, map(x->x*180/pi, theta_vals), af_db1, linewidth=2, label = "Progresssive Phase")
lines!(ax, map(x->x*180/pi, theta_vals), af_db2, linewidth=2, label = "Conjugate Symmetric")
ylims!(-80, 5)
for reg in p.sidelobe_regions
    points = 180/pi * map(x->x.θ, reg.region.points)
    values = ones(size(points))
    values .*= reg.upper
    values .= 20*log10.(values)
    lines!(ax, points, values, color = :black, linestyle = :dash, linewidth = 1.4)
end
axislegend(ax)
fig