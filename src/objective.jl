"""
Abstract supertype for all synthesis objectives.
"""
abstract type AbstractObjective end

"""
Abstract supertype for objectives solved with a single optimization model.
"""
abstract type DirectObjective <: AbstractObjective end

"""
    Feasible()

Find any excitation satisfying the pattern constraints.
"""
struct Feasible <: DirectObjective end

"""
    MaxAF([dirs])

Maximize the real array-factor response in the supplied direction or directions.
If no direction is provided, beam directions from the pattern are used.
"""
struct MaxAF{P} <: DirectObjective
    directions::Vector{P}
end

MaxAF() = MaxAF(AbstractDirection[])
MaxAF(dir) = MaxAF([direction(dir)])
MaxAF(dirs::AbstractVector) = MaxAF([direction(d) for d in dirs])

"""
    MinPower()

Minimize total excitation power.
"""
struct MinPower <: DirectObjective end

"""
    MinSLL(region; lower_bound = 0.0, upper_bound = -20.0dB)
    MinSLL(regions; lower_bound = 0.0, upper_bound = -20.0dB)

Minimize the sidelobe level over one or more regions.
"""
struct MinSLL{R, L} <: DirectObjective
    regions::Vector{R}
    lower_bound::L
    upper_bound::L
end

MinSLL(region::Region; lower_bound = 0.0, upper_bound = -20.0dB) = MinSLL([region]; lower_bound, upper_bound)
MinSLL(regions::AbstractVector{<:Region}; lower_bound = 0.0, upper_bound = -20.0dB) = MinSLL(collect(regions), lower_bound, upper_bound)

"""
    MinIntegratedPower(region)
    MinIntegratedPower(regions)

Minimize integrated array-factor power over one or more sampled regions.
"""
struct MinIntegratedPower{R} <: DirectObjective
    regions::Vector{R}
end

MinIntegratedPower(region::Region) = MinIntegratedPower([region])
MinIntegratedPower(regions::AbstractVector{<:Region}) = MinIntegratedPower{eltype(regions)}(collect(regions))


"""
    MinL1(; sum_limit = nothing)

Minimize a convex L1-type bound on the excitation variables.
"""
struct MinL1{L} <: DirectObjective
    sum_limit::L
end

MinL1(; sum_limit = nothing) = MinL1(sum_limit)

"""
    MinWeightedL1(alpha; sum_limit = nothing)

Minimize a weighted convex L1-type bound with one weight per excitation variable.
"""
struct MinWeightedL1{A, L} <: DirectObjective
    alpha::A
    sum_limit::L
end

MinWeightedL1(alpha; sum_limit = nothing) = MinWeightedL1(collect(alpha), sum_limit)


"""
    MinFieldError(region, reference)

Minimize squared array-factor error against a scalar or vector reference over `region`.
"""
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

# Wrapper for array_factor when used in model building context (with JuMP variables)
array_factor(model, array, points, weights, vars) = array_factor_reim(model, array, points, weights, vars)

objective!(model, objective, pattern, array, weights, vars, formulation) = error("Objective $(typeof(objective)) is not implemented.")

function objective!(model, ::Feasible, pattern, array, weights, vars, formulation)
    @objective(model, Min, 0.0)
    return nothing
end

function af_directions(objective::MaxAF, pattern)
    isempty(objective.directions) && return [b.direction for b in pattern.beams]
    return objective.directions
end

function objective!(model, objective::MaxAF, pattern, array, weights, vars, formulation)
    dirs = af_directions(objective, pattern)
    isempty(dirs) && error("MaxAF needs at least one direction.")
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

function region_sll_constraints!(model, region::Region, upper, array, weights, vars, formulation, robustness = nothing)
    af_re, af_im = array_factor(model, array, region.points, weights, vars)
    for i in eachindex(region.points)
        bound = robust_bound(upper, robustness, array, region.points[i])
        modulus_upper_bound!(model, af_re[i], imag_part(af_im, i), bound, formulation)
    end
end

function constrain_sll_objective!(model, objective::MinSLL, pattern, array, weights, vars, formulation, robustness = nothing)
    t = @variable(model, [1:length(objective.regions)])
    for (i, region) in enumerate(objective.regions)
        objective.lower_bound !== nothing && @constraint(model, t[i] >= objective.lower_bound)
        objective.lower_bound !== nothing && @constraint(model, t[i] <= objective.upper_bound)
        region_sll_constraints!(model, region, t[i], array, weights, vars, formulation, robustness)
    end
    return t
end

function objective!(model, objective::MinSLL, pattern, array, weights, vars, formulation)
    t = constrain_sll_objective!(model, objective, pattern, array, weights, vars, formulation)
    @objective(model, Min, sum(t))
    return t
end

function objective!(model, objective::MinSLL, pattern, array, weights, vars, formulation, robustness)
    t = constrain_sll_objective!(model, objective, pattern, array, weights, vars, formulation, robustness)
    @objective(model, Min, sum(t))
    return t
end

function objective!(model, objective, pattern, array, weights, vars, formulation, robustness)
    return objective!(model, objective, pattern, array, weights, vars, formulation)
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

array_factor_error_expression(re, im::Nothing, reference) = (re - real(reference))^2
array_factor_error_expression(re, im, reference) = (re - real(reference))^2 + (im - imag(reference))^2

function objective!(model, objective::MinFieldError, pattern, array, weights, vars, formulation)
    af_re, af_im = array_factor(model, array, objective.points, weights, vars)
    @objective(model, Min, sum(
        array_factor_error_expression(af_re[i], imag_part(af_im, i), objective.reference[i])
        for i in eachindex(objective.points)))
    return nothing
end

objective_beam_directions(obj::MaxAF, pattern) = af_directions(obj, pattern)
objective_beam_directions(::AbstractObjective, pattern) = AbstractDirection[]

check_formulation(obj::Union{MinPower, MinIntegratedPower, MinFieldError}, formulation) = formulation isa LP && error("$(typeof(obj)) has a quadratic objective, incompatible with LP.")
check_formulation(::AbstractObjective, formulation) = nothing
