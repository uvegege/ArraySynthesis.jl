abstract type AbstractObjective end
abstract type DirectObjective <: AbstractObjective end

struct Feasible <: DirectObjective end

struct MaxGain{P} <: DirectObjective
    directions::Vector{P}
end

MaxGain() = MaxGain(AbstractDirection[])
MaxGain(dir) = MaxGain([direction(dir)])
MaxGain(dirs::AbstractVector) = MaxGain([direction(d) for d in dirs])

struct MinPower <: DirectObjective end

struct MinSLL{R, L} <: DirectObjective
    regions::Vector{R}
    lower_bound::L
    upper_bound::L
end

MinSLL(region::Region; lower_bound = 0.0, upper_bound = -20.0dB) = MinSLL([region]; lower_bound, upper_bound)
MinSLL(regions::AbstractVector{<:Region}; lower_bound = 0.0, upper_bound = -20.0dB) = MinSLL(collect(regions), lower_bound, upper_bound)

struct MinIntegratedPower{R} <: DirectObjective
    regions::Vector{R}
end

MinIntegratedPower(region::Region) = MinIntegratedPower([region])
MinIntegratedPower(regions::AbstractVector{<:Region}) = MinIntegratedPower{eltype(regions)}(collect(regions))


struct MinL1{L} <: DirectObjective
    sum_limit::L
end

MinL1(; sum_limit = nothing) = MinL1(sum_limit)

struct MinWeightedL1{A, L} <: DirectObjective
    alpha::A
    sum_limit::L
end

MinWeightedL1(alpha; sum_limit = nothing) = MinWeightedL1(collect(alpha), sum_limit)


struct MinActiveElements <: DirectObjective end

struct MinFieldError{P, T} <: DirectObjective
    points::Vector{P}
    reference::Vector{T}
end

function MinFieldError(region::Region, reference::AbstractVector)
    length(reference) == length(region.points) || error("Reference length must match region points.")
    return MinFieldError(region.points, collect(reference))
end

function MinFieldError(region::Region, reference::Number)
    return MinFieldError(region.points, fill(reference, length(region.points)))
end

# Minimize integrated power subject to Re(AF(θ₀)) = 1, which maximizes
# directivity = |AF(θ₀)|² / ∫|AF(θ)|² dΩ. SOCP only.
struct MaxDirectivity{P, R} <: DirectObjective
    direction::P
    region::R
end

MaxDirectivity(dir, region::Region) = MaxDirectivity(direction(dir), region)

# Wrapper for array_factor when used in model building context (with JuMP variables)
array_factor(model, array, points, weights, vars) = array_factor_reim(model, array, points, weights, vars)

objective!(model, objective, pattern, array, weights, vars, formulation) = error("Objective $(typeof(objective)) is not implemented.")

function objective!(model, ::Feasible, pattern, array, weights, vars, formulation)
    @objective(model, Min, 0.0)
    return nothing
end

function gain_directions(objective::MaxGain, pattern)
    isempty(objective.directions) && return [b.direction for b in pattern.beams]
    return objective.directions
end

function objective!(model, objective::MaxGain, pattern, array, weights, vars, formulation)
    dirs = gain_directions(objective, pattern)
    isempty(dirs) && error("MaxGain needs at least one direction.")
    af_re, af_im = array_factor(model, array, dirs, weights, vars)
    @objective(model, Max, sum(af_re))
    return nothing
end

power_expression(vars::SparseVariables) = power_expression(vars.variables)
power_expression(vars::AmplitudeVariables) = sum(vars.a[n]^2 for n in eachindex(vars.a))

function power_expression(vars::WeightVariables)
    return sum(vars.w_re[n]^2 for n in eachindex(vars.w_re)) + sum(vars.w_im[n]^2 for n in eachindex(vars.w_im))
end

function objective!(model, ::MinPower, pattern, array, weights, vars, formulation)
    @objective(model, Min, power_expression(vars))
    return nothing
end

pattern_power_expression(re, im::Nothing) = sum(re[i]^2 for i in eachindex(re))
pattern_power_expression(re, im) = sum(re[i]^2 + im[i]^2 for i in eachindex(re))

function objective!(model, objective::MinIntegratedPower, pattern, array, weights, vars, formulation)
    power = 0.0
    for region in objective.regions
        af_re, af_im = array_factor(model, array, region.points, weights, vars)
        power += pattern_power_expression(af_re, af_im)
    end
    @objective(model, Min, power)
    return nothing
end

function region_sll_constraints!(model, region::Region, upper, array, weights, vars, formulation)
    af_re, af_im = array_factor(model, array, region.points, weights, vars)
    for i in eachindex(region.points)
        modulus_upper_bound!(model, af_re[i], imag_part(af_im, i), upper, formulation)
    end
end

function constrain_sll_objective!(model, objective::MinSLL, pattern, array, weights, vars, formulation)
    t = @variable(model, [1:length(objective.regions)])
    for (i, region) in enumerate(objective.regions)
        objective.lower_bound !== nothing && @constraint(model, t[i] >= objective.lower_bound)
        objective.lower_bound !== nothing && @constraint(model, t[i] <= objective.upper_bound)
        region_sll_constraints!(model, region, t[i], array, weights, vars, formulation)
    end
    return t
end

function objective!(model, objective::MinSLL, pattern, array, weights, vars, formulation)
    t = constrain_sll_objective!(model, objective, pattern, array, weights, vars, formulation)
    @objective(model, Min, sum(t))
    return t
end

objective_nvariables(vars::SparseVariables) = nvariables(vars.variables)
objective_nvariables(vars) = nvariables(vars)


function objective!(model, objective::MinL1, pattern, array, weights, vars, formulation)
    t = l1_variables!(model, vars)
    l1_bound!(model, t, vars, formulation)
    l1_limit!(model, t, objective.sum_limit)
    @objective(model, Min, sum(t))
    return t
end

function objective!(model, objective::MinWeightedL1, pattern, array, weights, vars, formulation)
    objective_nvariables(vars) == length(objective.alpha) || error("Alpha length must match the number of variables.")
    t = l1_variables!(model, vars)
    l1_bound!(model, t, vars, formulation)
    l1_limit!(model, t, objective.sum_limit)
    @objective(model, Min, sum(objective.alpha[n] * t[n] for n in eachindex(t)))
    return t
end

#=
function objective!(model, ::MinActiveElements, pattern, array, weights, vars::SparseVariables, formulation::MILP)
    @objective(model, Min, sum(vars.active))
    return vars.active
end

function objective!(model, ::MinActiveElements, pattern, array, weights, vars, formulation)
    error("MinActiveElements requires MILP variables.")
end
=#

field_error_expression(re, im::Nothing, reference) = (re - real(reference))^2
field_error_expression(re, im, reference) = (re - real(reference))^2 + (im - imag(reference))^2

function objective!(model, objective::MinFieldError, pattern, array, weights, vars, formulation)
    af_re, af_im = array_factor(model, array, objective.points, weights, vars)
    @objective(model, Min, sum(
        field_error_expression(af_re[i], imag_part(af_im, i), objective.reference[i])
        for i in eachindex(objective.points)))
    return nothing
end

function objective!(model, objective::MaxDirectivity, pattern, array, weights, vars, formulation::SOCP)
    af_re, af_im = array_factor(model, array, [objective.direction], weights, vars)
    @constraint(model, af_re[1] == 1.0)
    paf_re, paf_im = array_factor(model, array, objective.region.points, weights, vars)
    @objective(model, Min, pattern_power_expression(paf_re, paf_im))
    return nothing
end

function objective!(model, ::MaxDirectivity, pattern, array, weights, vars, formulation)
    error("MaxDirectivity requires SOCP formulation.")
end

objective_beam_directions(obj::MaxGain, pattern) = gain_directions(obj, pattern)
objective_beam_directions(::AbstractObjective, pattern) = AbstractDirection[]

check_formulation(::MaxDirectivity, formulation) = formulation isa SOCP || error("MaxDirectivity requires SOCP.")
#check_formulation(::MinActiveElements, formulation) = formulation isa MILP || error("MinActiveElements requires MILP.")
check_formulation(obj::Union{MinPower, MinIntegratedPower, MinFieldError}, formulation) = formulation isa LP && error("$(typeof(obj)) has a quadratic objective, incompatible with LP.")
check_formulation(::AbstractObjective, formulation) = nothing
