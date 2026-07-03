abstract type AbstractExcitationModel end
abstract type AbstractFormulation end

struct LP <: AbstractFormulation
    polygon_faces::Int
end
LP(; polygon_faces::Int = 8) = LP(polygon_faces)

struct QP <: AbstractFormulation
    polygon_faces::Int
end
QP(; polygon_faces::Int = 8) = QP(polygon_faces)

struct SOCP <: AbstractFormulation end

struct MILP <: AbstractFormulation
    max_active_elements::Union{Nothing, Int}
    big_m::Float64
    polygon_faces::Int
end

function MILP(; max_active_elements = nothing, big_m::Real = 10.0, polygon_faces::Int = 8)
    MILP(max_active_elements, Float64(big_m), polygon_faces)
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

#=
function variables!(model, array, weights, formulation::MILP)
    variables = variables!(model, array, weights)
    active = @variable(model, [1:nvariables(variables)], Bin)
    bind_activity!(model, variables, active, formulation.big_m)
    formulation.max_active_elements !== nothing && @constraint(model, sum(active) <= formulation.max_active_elements)
    return SparseVariables(variables, active)
end

function bind_activity!(model, x::WeightVariables, active, big_m)
    for n in eachindex(active)
        @constraint(model, x.w_re[n] <= big_m * active[n])
        @constraint(model, x.w_re[n] >= -big_m * active[n])
    end
    for n in eachindex(x.w_im)
        @constraint(model, x.w_im[n] <= big_m * active[n])
        @constraint(model, x.w_im[n] >= -big_m * active[n])
    end
end

function bind_activity!(model, x::AmplitudeVariables, active, big_m)
    for n in eachindex(active)
        @constraint(model, x.a[n] <= big_m * active[n])
        @constraint(model, x.a[n] >= -big_m * active[n])
    end
end
=#
