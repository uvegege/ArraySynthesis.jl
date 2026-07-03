using LinearAlgebra: factorize

direction_vector(p::ThetaDirection) =(sin(p.θ), 0.0, cos(p.θ))
direction_vector(p::UVDirection) = (p.u, p.v, sqrt(max(0.0, 1.0 - p.u^2 - p.v^2)))

function direction_matrix(points)
    U = zeros(Float64, 3, length(points))
    for (i, p) in enumerate(points)
        U[:, i] .= direction_vector(p)
    end
    return U
end

array_factor_reim(model, array, points, weights, vars::SparseVariables) = array_factor_reim(model, array, points, weights, vars.variables)

imag_part(im::Nothing, i) = nothing
imag_part(im, i) = im[i]

steering_matrix_for(array, points, weights) = steering_matrix(array, direction_matrix(points))
steering_matrix_for(array, points, weights::ProgressivePhaseAmplitude) = steering_matrix(array, direction_matrix(points) .- phase_direction(weights.β))

function array_factor_reim_from_steering(A, ::ComplexWeights, vars)
    A_re = real.(A); A_im = imag.(A)
    af_re = A_re * vars.w_re - A_im * vars.w_im
    af_im = A_im * vars.w_re + A_re * vars.w_im
    return af_re, af_im
end

function array_factor_reim(model, array, points, weights::ComplexWeights, vars)
    A = steering_matrix_for(array, points, weights)
    return array_factor_reim_from_steering(A, weights, vars)
end

function array_factor_reim_from_steering(A, ::RealWeights, vars::WeightVariables)
    af_re = real.(A) * vars.w_re
    af_im = imag.(A) * vars.w_re
    return af_re, af_im
end

function array_factor_reim(model, array, points, weights::RealWeights, vars::WeightVariables)
    A = steering_matrix_for(array, points, weights)
    return array_factor_reim_from_steering(A, weights, vars)
end

function array_factor_reim_from_steering(A::Tuple, ::ConjugateSymmetricWeights, vars::WeightVariables)
    A_cos, A_sin = A
    af_re = A_cos * vars.w_re + A_sin * vars.w_im
    return af_re, nothing
end

function array_factor_reim(model, array::SymmetricArray, points, weights::ConjugateSymmetricWeights, vars::WeightVariables)
    A = steering_matrix_for(array, points, weights)
    return array_factor_reim_from_steering(A, weights, vars)
end

function array_factor_reim_from_steering(A, ::ProgressivePhaseAmplitude, vars::AmplitudeVariables)
    af_re = real.(A) * vars.a
    af_im = imag.(A) * vars.a
    return af_re, af_im
end

function array_factor_reim(model, array::ArrayGeometry, points, weights::ProgressivePhaseAmplitude, vars::AmplitudeVariables)
    A = steering_matrix_for(array, points, weights)
    return array_factor_reim_from_steering(A, weights, vars)
end

function array_factor_reim_from_steering(A::Tuple, ::ProgressivePhaseAmplitude, vars::AmplitudeVariables)
    A_cos, _ = A
    return A_cos * vars.a, nothing
end

function array_factor_reim(model, array::SymmetricArray, points, weights::ProgressivePhaseAmplitude, vars::AmplitudeVariables)
    A = steering_matrix_for(array, points, weights)
    return array_factor_reim_from_steering(A, weights, vars)
end


phase_direction(p::ThetaDirection) = (sin(p.θ), 0.0, cos(p.θ))
phase_direction(β::UVDirection) = (β.u, β.v, sqrt(1.0 - β.u^2 - β.v^2))
phase_direction(β::Number) = (sin(β), 0.0, cos(β))
phase_direction(β::Tuple{<:Real, <:Real}) = (β[1], β[2], sqrt(1.0 - β[1]^2 - β[2]^2))

#phase_direction(β::AbstractVector{<:Real}) = length(β) == 3 ? collect(β) : error("β vector must have length 3.")

function polygon_bound!(model, re, im, bound, faces)
    for l in 0:(faces - 1)
        θ = 2π * l / faces
        @constraint(model, cos(θ) * re - sin(θ) * im <= bound * cos(π / faces))
    end
end

function polygon_bound!(model, re::AffExpr, im::AffExpr, bound, faces)
    rhs = bound * cos(π / faces)
    for l in 0:(faces - 1)
        θ = 2π * l / faces
        c = cos(θ); s = sin(θ)
        lhs = AffExpr(c * re.constant - s * im.constant)
        for (var, coef) in re.terms
            add_to_expression!(lhs, c * coef, var)
        end
        for (var, coef) in im.terms
            add_to_expression!(lhs, -s * coef, var)
        end
        @constraint(model, lhs <= rhs)
    end
end

bound_at(bound::AbstractVector, i) = bound[i]
bound_at(bound, i) = bound

function modulus_upper_bound!(model, re, im::Nothing, bound, formulation::Union{LP, QP, MILP})
    @constraint(model, re <= bound)
    @constraint(model, re >= -bound)
end

function modulus_upper_bound!(model, re, im::Nothing, bound, formulation::SOCP)
    @constraint(model, re <= bound)
    @constraint(model, re >= -bound)
end

function modulus_upper_bound!(model, re, im, bound, formulation::Union{LP, QP, MILP})
    polygon_bound!(model, re, im, bound, formulation.polygon_faces)
end

function modulus_upper_bound!(model, re, im, bound, formulation::SOCP)
    @constraint(model, [bound, re, im] in SecondOrderCone())
end

function l1_variables!(model, vars)
    return @variable(model, [1:nvariables(vars)])
end

function l1_variables!(model, vars::SparseVariables)
    return l1_variables!(model, vars.variables)
end

function l1_bound!(model, t, vars::SparseVariables, formulation)
    l1_bound!(model, t, vars.variables, formulation)
end

function l1_bound!(model, t, vars::WeightVariables, formulation::Union{LP, QP, MILP})
    for n in eachindex(t)
        @constraint(model, t[n] >= vars.w_re[n])
        @constraint(model, t[n] >= -vars.w_re[n])
    end
    for n in eachindex(vars.w_im)
        @constraint(model, t[n] >= vars.w_im[n])
        @constraint(model, t[n] >= -vars.w_im[n])
    end
end

function l1_bound!(model, t, vars::WeightVariables, formulation::SOCP)
    for n in eachindex(vars.w_im)
        @constraint(model, [t[n], vars.w_re[n], vars.w_im[n]] in SecondOrderCone())
    end
    for n in (length(vars.w_im) + 1):length(t)
        @constraint(model, t[n] >= vars.w_re[n])
        @constraint(model, t[n] >= -vars.w_re[n])
    end
end

function l1_bound!(model, t, vars::AmplitudeVariables, formulation)
    for n in eachindex(t)
        @constraint(model, t[n] >= vars.a[n])
    end
end

function shared_l1_bound!(model, t, variables, formulation)
    for vars in variables
        l1_bound!(model, t, vars, formulation)
    end
end

#=
function active_big_m!(model, vars::SparseVariables, big_m)
    active_big_m!(model, vars.variables, vars.active, big_m)
end

function active_big_m!(model, vars::WeightVariables, active, big_m)
    for n in eachindex(active)
        @constraint(model, vars.w_re[n] <= big_m * active[n])
        @constraint(model, vars.w_re[n] >= -big_m * active[n])
    end
    for n in eachindex(vars.w_im)
        @constraint(model, vars.w_im[n] <= big_m * active[n])
        @constraint(model, vars.w_im[n] >= -big_m * active[n])
    end
end

function active_big_m!(model, vars::AmplitudeVariables, active, big_m)
    for n in eachindex(active)
        @constraint(model, vars.a[n] <= big_m * active[n])
    end
end

active_limit!(model, vars::SparseVariables, limit::Nothing) = nothing
active_limit!(model, vars::SparseVariables, limit::Integer) = @constraint(model, sum(vars.active) <= limit)
=#

l1_limit!(model, t, limit::Nothing) = nothing
l1_limit!(model, t, limit) = @constraint(model, sum(t) <= limit)


function sidelobe_constraints!(model, sl::SideLobeRegion, upper, array, weights, vars, formulation)
    af_re, af_im = array_factor_reim(model, array, sl.region.points, weights, vars)
    for i in eachindex(sl.region.points)
        modulus_upper_bound!(model, af_re[i], imag_part(af_im, i), bound_at(upper, i), formulation)
    end
end

function sidelobe_constraints!(model, sl::SideLobeRegion, array, weights, vars, formulation)
    sidelobe_constraints!(model, sl, sl.upper, array, weights, vars, formulation)
end

shaped_reference(array, points, weights, target, im::Nothing) = target
shaped_reference(array, points, weights, target, im) = phase_retrieve(reference_matrix(array, points, weights), target)
shaped_reference_from_steering(A, target, im::Nothing) = target
shaped_reference_from_steering(A, target, im) = phase_retrieve(A, target)


reference_matrix(array, points, weights) = steering_matrix(array, direction_matrix(points))

reference_matrix(array, points, weights::ProgressivePhaseAmplitude) = steering_matrix(array, direction_matrix(points) .- phase_direction(weights.β))


function phase_retrieve(A, target; iterations = 100)
    y = ComplexF64.(target, 0.0)
    F = factorize(A)
    for _ in 1:iterations
        w = F \ y
        y .= abs.(y) .* cis.(angle.(A * w))
    end
    return y
end

function shaped_bound!(model, re, im::Nothing, target, ripple, formulation::Union{LP, QP, MILP})
    b1 = target * ripple
    b2 = target / ripple
    @constraint(model, re >= min(b1, b2))
    @constraint(model, re <= max(b1, b2))
end

function shaped_bound!(model, re, im::Nothing, target, ripple, ::SOCP)
    b1 = target * ripple
    b2 = target / ripple
    @constraint(model, re >= min(b1, b2))
    @constraint(model, re <= max(b1, b2))
end

function shaped_bound!(model, re, im, target, ripple, formulation::Union{LP, QP, MILP})
    ripple_radius = abs(target) * abs(1 - ripple)
    polygon_bound!(model, re - real(target), im - imag(target), ripple_radius, formulation.polygon_faces)
end

function shaped_bound!(model, re, im, target, ripple, ::SOCP)
    ε = abs(target) * abs(1 - ripple)
    @constraint(model, [ε, re - real(target), im - imag(target)] in SecondOrderCone())
end

function gain_constraint!(model, re, im::Nothing, gain)
    @constraint(model, re == gain)
end

function gain_constraint!(model, re, im, gain)
    @constraint(model, re == gain)
    @constraint(model, im == 0)
end

constrain_variables!(model, array, weights, vars) = nothing
constrain_variables!(model, array, weights, vars::SparseVariables) = constrain_variables!(model, array, weights, vars.variables)

function constrain_variables!(model, array::SymmetricArray, ::ConjugateSymmetricWeights, vars::WeightVariables)
    for n in axes(array.positions, 2)
        is_origin(view(array.positions, :, n)) && @constraint(model, vars.w_im[n] == 0)
    end
end

function constraints!(model, pattern::Pattern, array, weights, vars, formulation)
    constrain_variables!(model, array, weights, vars)

    for b in pattern.beams
        af_re, af_im = array_factor_reim(model, array, (b.direction, ), weights, vars)
        gain_constraint!(model, af_re[1], imag_part(af_im, 1), b.gain)
    end

    for sb in pattern.shaped_beams
        A = steering_matrix_for(array, sb.region.points, weights)
        af_re, af_im = array_factor_reim_from_steering(A, weights, vars)
        target = shaped_reference_from_steering(A, sb.target, af_im)
        for i in eachindex(sb.region.points)
            shaped_bound!(model, af_re[i], imag_part(af_im, i), target[i], sb.ripple, formulation)
        end
    end

    for n in pattern.null_directions
        af_re, af_im = array_factor_reim(model, array, (n.direction, ), weights, vars)
        modulus_upper_bound!(model, af_re[1], imag_part(af_im, 1), n.level, formulation)
    end

    for sl in pattern.sidelobe_regions
        sidelobe_constraints!(model, sl, array, weights, vars, formulation)
    end
end
