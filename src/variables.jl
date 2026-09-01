abstract type AbstractExcitationModel end
"""
Abstract supertype for optimization formulations.
"""
abstract type AbstractFormulation end

"""
    LP(; polygon_faces = 8)

Linear-programming formulation. Complex magnitude bounds are approximated by a
regular polygon with `polygon_faces` sides.
"""
struct LP <: AbstractFormulation
    polygon_faces::Int
end
LP(; polygon_faces::Int = 8) = LP(polygon_faces)

"""
    QP(; polygon_faces = 8)

Quadratic-programming formulation with polygonal magnitude bounds.
"""
struct QP <: AbstractFormulation
    polygon_faces::Int
end
QP(; polygon_faces::Int = 8) = QP(polygon_faces)

"""
    SOCP()

Second-order cone formulation using exact conic magnitude bounds.
"""
struct SOCP <: AbstractFormulation end

"""
    MILP(; big_m = 10.0, polygon_faces = 8)

Mixed-integer formulation for sparse active-element constraints.
"""
struct MILP <: AbstractFormulation
    big_m::Union{Float64, Vector{Float64}}
    polygon_faces::Int
end

function MILP(;big_m = 10.0, polygon_faces::Int = 8)
    MILP(Float64.(big_m), polygon_faces)
end

struct WeightVariables{R,I}
    w_re::R
    w_im::I
end

struct AmplitudeVariables{A}
    a::A
end

nvariables(x::WeightVariables) = length(x.w_re)
nvariables(x::AmplitudeVariables) = length(x.a)

struct SparseVariables{V, B}
    variables::V
    active::B
end 

function variables!(model, array, weights, formulation::Union{LP, QP, SOCP})
    return variables!(model, array, weights)
end

function variables!(model, array::ArrayGeometry, ::RealWeights)
    N = size(array.positions, 2)
    w_re = @variable(model, [1:N])
    w_im = w_re[1:0]
    return WeightVariables(w_re, w_im)
end


function variables!(model, array::SymmetricArray, ::ComplexWeights)
    error("These weights do not use symmetric representatives. Use materialize(array) explicitly.")
end

function variables!(model, array::SymmetricArray, ::RealWeights)
    error("These weights do not use symmetric representatives. Use materialize(array) explicitly.")
end

function variables!(model, ::ArrayGeometry, ::ConjugateSymmetricWeights)
    error("ConjugateSymmetricWeights requires a SymmetricArray.")
end

function variables!(model, array, ::Union{ComplexWeights, ConjugateSymmetricWeights})
    N = size(array.positions, 2)
    w_re = @variable(model, [1:N])
    w_im = @variable(model, [1:N])
    return WeightVariables(w_re, w_im)
end

real_weight(vars::WeightVariables, n) = vars.w_re[n]
imag_weight(vars::WeightVariables, n) = vars.w_im[n]

function variables!(model, array, ::ProgressivePhaseAmplitude)
    N = size(array.positions, 2)
    #a = @variable(model, [1:N], lower_bound = 0)
    a = @variable(model, [1:N])
    return AmplitudeVariables(a)
end

variables!(model, array, ::QuantizedAmplitude, formulation::Union{LP, QP, SOCP}) = error("QuantizedAmplitude requires the MILP formulation.")
variables!(model, array, ::QuantizedPhase, formulation::Union{LP,QP,SOCP}) = error("QuantizedPhase requires the MILP formulation.")

function variables!(model, array, weights::QuantizedAmplitude, formulation::MILP)
    N = size(array.positions, 2)

    levels = weights.levels
    M = length(levels)
    U = formulation.big_m
    U isa Float64 ||  error("QuantizedAmplitude needs a global (scalar) big_m.")
    
    M >= 2 || error("QuantizedAmplitude needs at least two levels.")
    #all(>=(0), levels) || error("QuantizedAmplitude levels must be nonnegative.") #TODO CHECK

    lvls = sort(collect(Float64, levels))
    Δ = lvls[2] - lvls[1]
    are_uniform = all(k -> isapprox(lvls[k], lvls[1] + (k - 1) * Δ; atol = 1e-9), eachindex(lvls))
    lmin = lvls[1]
    
    if are_uniform && !weights.relative
        K = length(lvls) - 1
        q = @variable(model, [1:N], integer = true, lower_bound = 0, upper_bound = K)
        return AmplitudeVariables([lmin + Δ * q[n] for n in 1:N])
    end

    if are_uniform && weights.relative
        isapprox(maximum(abs, levels), 1.0) || error("Relative QuantizedAmplitude levels must be normalized so their maximum absolute value is 1.0.")
        K = length(lvls) - 1
        #U = formulation.big_m # TODO: No estoy seguro de si este valor es el más razonable
        n_B = ceil(Int, log2(K + 1))
        powers = [2^(b - 1) for b in 1:n_B]
        y = @variable(model, [1:N, 1:n_B], Bin)
        @constraint(model, [n in 1:N], sum(powers[b] * y[n, b] for b in 1:n_B) <= K)
        V = @variable(model, lower_bound = 0, upper_bound = U)
        s = @variable(model, [1:N, 1:n_B], lower_bound = 0)
        @constraint(model, [n in 1:N, b in 1:n_B], s[n, b] <= U * y[n, b])
        @constraint(model, [n in 1:N, b in 1:n_B], s[n, b] <= V)
        @constraint(model, [n in 1:N, b in 1:n_B], s[n, b] >= V - U * (1 - y[n, b]))
        return AmplitudeVariables([lmin * V + Δ * sum(powers[b] * s[n, b] for b in 1:n_B) for n in 1:N])
    end

    z = @variable(model, [1:N, 1:M], Bin)
    @constraint(model, [n in 1:N], sum(z[n, k] for k in 1:M) == 1)

    if !weights.relative
        return AmplitudeVariables([sum(levels[k] * z[n, k] for k in 1:M) for n in 1:N])
    end

    isapprox(maximum(abs, levels), 1.0) || error("Relative QuantizedAmplitude levels must be normalized so their maximum is 1.0.")

    V = @variable(model, lower_bound = 0) 
    v = @variable(model, [1:N, 1:M], lower_bound = 0)

    @constraint(model, [n in 1:N, k in 1:M], v[n, k] <= U * z[n, k])
    @constraint(model, [n in 1:N], sum(v[n, k] for k in 1:M) == V)

    return AmplitudeVariables([sum(levels[k] * v[n, k] for k in 1:M) for n in 1:N])
end

function variables!(model, array, weights::QuantizedPhase, formulation::MILP)
    N = size(array.positions, 2)
    discrete_phases = weights.levels
    M = length(discrete_phases)
    M >= 2 || error("QuantizedPhase needs at least two phases.")
    c = cos.(discrete_phases)
    s = sin.(discrete_phases)

    z = @variable(model, [1:N, 1:M], Bin)
    @constraint(model, [n in 1:N], sum(z[n, k] for k in 1:M) == 1)

    amp = weights.amplitude
    amp_n(n) = amp isa AbstractVector ? amp[n] : amp

    if !weights.free_amplitudes
        w_re = [amp_n(n) * sum(z[n, k] * c[k] for k in 1:M) for n in 1:N]
        w_im = [amp_n(n) * sum(z[n, k] * s[k] for k in 1:M) for n in 1:N]
        return WeightVariables(w_re, w_im)
    end

    U = formulation.big_m
    U_n(n) = U isa AbstractVector ? U[n] : U

    v = @variable(model, [1:N, 1:M], lower_bound = 0)
    @constraint(model, [n in 1:N, k in 1:M], v[n, k] <= U_n(n) * z[n, k])
    w_re = [sum(v[n, k] * c[k] for k in 1:M) for n in 1:N]
    w_im = [sum(v[n, k] * s[k] for k in 1:M) for n in 1:N]
    return WeightVariables(w_re, w_im)
end

