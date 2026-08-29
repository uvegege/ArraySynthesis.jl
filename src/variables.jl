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
    big_m::Float64
    polygon_faces::Int
end

function MILP(;big_m::Real = 100.0, polygon_faces::Int = 8)
    MILP(Float64(big_m), polygon_faces)
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

variables!(model, array, ::QuantizedAmplitude, formulation::Union{LP, QP, SOCP}) = error("QuantizedAmplitude requires the MILP/MISOCP formulation.")

function variables!(model, array, weights::QuantizedAmplitude, formulation::MILP)
    N = size(array.positions, 2)

    levels = weights.levels
    M = length(levels)

    M >= 1 || error("QuantizedAmplitude needs at least one level.")
    #all(>=(0), levels) || error("QuantizedAmplitude levels must be nonnegative.") #TODO CHECK

    # Check if levels are uniform. 
    lvls = sort(collect(Float64, levels))
    Δ = lvls[2] - lvls[1]
    are_uniform = all(k -> isapprox(lvls[k], lvls[1] + (k - 1) * Δ; atol = 1e-9), eachindex(lvls))

    if are_uniform && !weights.relative
        lmin = lvls[1]
        K = length(lvls) - 1
        q = @variable(model, [1:N], integer = true, lower_bound = 0, upper_bound = K)
        return AmplitudeVariables([lmin + Δ * q[n] for n in 1:N])
    end

    z = @variable(model, [1:N, 1:M], Bin)
    @constraint(model, [n in 1:N], sum(z[n, k] for k in 1:M) == 1)

    if !weights.relative
        return AmplitudeVariables([sum(levels[k] * z[n, k] for k in 1:M) for n in 1:N])
    end

    isapprox(maximum(abs, levels), 1.0) || error("Relative QuantizedAmplitude levels must be normalized so their maximum is 1.0.")

    big_m = formulation.big_m
    V = @variable(model, lower_bound = 0) #V = @variable(model, 0 <= V <= big_m)
    v = @variable(model, [1:N, 1:M], lower_bound = 0) #v = @variable(model, [1:N, 1:M] >= 0)

    @constraint(model, [n in 1:N, k in 1:M], v[n, k] <= big_m * z[n, k])
    @constraint(model, [n in 1:N], sum(v[n, k] for k in 1:M) == V)

    return AmplitudeVariables([sum(levels[k] * v[n, k] for k in 1:M) for n in 1:N])
end
