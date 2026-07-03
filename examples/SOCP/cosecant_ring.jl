using ArraySynthesis
using ArraySynthesis: °, dB
using GLMakie
using LinearAlgebra
using Clarabel

array = circular_array(
    [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5],
    [1, 6, 12, 18, 25, 31, 37, 43, 50, 56],
)

function cosecbeam(ui::Real, vmax::Real, vmin::Real, amax::Real, amin::Real, ancho_umax::Real, ancho_umin::Real)
    n = 10
    umax_lin = collect(range(ui - ancho_umax / 2, ui + ancho_umax / 2, length = n))
    umin_lin = collect(range(ui - ancho_umin / 2, ui + ancho_umin / 2, length = n))
    ypos = collect(range(vmax, vmin, length = n))

    posu = Matrix{Float64}(undef, n, n)
    posv = Matrix{Float64}(undef, n, n)
    amp = Matrix{Float64}(undef, n, n)

    for v_idx in 1:n
        mu = (umin_lin[v_idx] - umax_lin[v_idx]) / (vmin - vmax)
        for u_idx in 1:n
            posu[u_idx, v_idx] = mu * (ypos[u_idx] - vmax) + umax_lin[v_idx]
            posv[u_idx, v_idx] = ypos[u_idx]
            amp[u_idx, v_idx] = csc(0.1 + 0.05 * u_idx)
        end
    end

    amp ./= maximum(amp)
    return vec(posu), vec(posv), vec(amp), vec(0.8 .* amp)
end

u_haz, v_haz, hsup, hinf = cosecbeam(0.0, 0.25, -0.1, 1.0, 0.2, 0.25, 0.45)
beam_points = [UVDirection(u_haz[i], v_haz[i]) for i in eachindex(u_haz)]
beam_region = Region(beam_points, :cosecant_beam)

cosecant_target = hsup

beam_shape = Polygon(
    [-0.125, 0.125, 0.225, -0.225],
    [0.25, 0.25, -0.1, -0.1],
)
bp = 0.2
uv_step = 4°
sll_region = visible_region(beam_shape; step = uv_step, bandpass = bp, filtered = false)

p = pattern(
    shaped_beam(beam_region, cosecant_target, ripple = 0.8),
    sidelobes(sll_region, -22dB),
)

coef = ComplexWeights()
obj = MinSLL([sll_region])
result = synthesize(array, p, obj, coef, SOCP(), Clarabel.Optimizer)

U = collect(-1.0:0.02:1.0)
V = collect(-1.0:0.02:1.0)
dirs = [UVDirection(u, v) for u in U for v in V]
AF = array_factor(array, coef, result.weights, dirs)
af_db = 20 .* log10.(max.(abs.(AF), 1e-12))
af_vals = reshape(af_db, length(V), length(U))'

begin
    fig = Figure()
    U = collect(-1.0:0.01:1.0)
    V = collect(-1.0:0.01:1.0)
    dirs = [UVDirection(u,v) for u in U, v in V]
    AF = array_factor(array, coef, result.weights, dirs)
    af_vals = reshape(20 .* log10.(abs.(getindex.(AF, 1)) .+ 1e-12), length(U), length(V))
    af_vals = [max(af_i, -50) for af_i in af_vals]
    ax3 = Axis3(fig[1,1])
    surface!(ax3, U, V, af_vals, colormap = parula_cm)
    fig
end
