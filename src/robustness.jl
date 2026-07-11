struct Tolerances{P}
    pointing::P
    phase::Float64
    amplitude::Float64
    position::NTuple{3, Float64}
    combine::Symbol
end

struct RobustMargin{T, R}
    tolerances::T
    norm_bound::R
end

function robust(; pointing_accuracy = 0.0, phase_tolerance = 0.0, amplitude_tolerance = 0.0, position_tolerance = 0.0, combine = :quadrature)
    combine in (:quadrature, :sum) || throw(ArgumentError("combine must be :quadrature or :sum."))
    return Tolerances(pointing_tuple(pointing_accuracy), Float64(phase_tolerance), amplitude_linear_tolerance(amplitude_tolerance), position_tuple(position_tolerance), combine)
end

pointing_tuple(x::Number) = (Float64(x), Float64(x))
pointing_tuple(x::Tuple{<:Real, <:Real}) = (Float64(x[1]), Float64(x[2]))
pointing_tuple(x::NamedTuple{(:u, :v)}) = (Float64(x.u), Float64(x.v))

function amplitude_linear_tolerance(x::Real)
    y = Float64(x)
    isapprox(y, 1.0; atol = 10eps(Float64), rtol = 0.0) && return 0.0
    return y > 0.5 ? abs(y - 1) : y
end

position_tuple(x::Number) = (Float64(x), Float64(x), Float64(x))
position_tuple(x::Tuple{<:Real, <:Real, <:Real}) = (Float64(x[1]), Float64(x[2]), Float64(x[3]))
position_tuple(x::NamedTuple{(:xy, :z)}) = (Float64(x.xy), Float64(x.xy), Float64(x.z))
position_tuple(x::NamedTuple{(:x, :y, :z)}) = (Float64(x.x), Float64(x.y), Float64(x.z))

direction_vector(p::ThetaDirection) = (sin(p.θ), 0.0, cos(p.θ))
direction_vector(p::UVDirection) = (p.u, p.v, sqrt(max(0.0, 1.0 - p.u^2 - p.v^2)))

theta_direction_derivative(p::ThetaDirection) = (cos(p.θ), 0.0, -sin(p.θ))

function uv_direction_derivatives(p::UVDirection)
    w = sqrt(max(eps(Float64), 1.0 - p.u^2 - p.v^2))
    return (1.0, 0.0, -p.u / w), (0.0, 1.0, -p.v / w)
end

physical_array(array::ArrayGeometry) = array
physical_array(array::SymmetricArray) = materialize(array)

function steering_derivative_norm(array, point::ThetaDirection)
    du = theta_direction_derivative(point)
    total = 0.0
    for p in eachcol(array.positions)
        coeff = 2*π * dot(p, du)
        total += coeff^2
    end
    return sqrt(total)
end

function steering_derivative_matrix(array, point::UVDirection)
    du, dv = uv_direction_derivatives(point)
    J = zeros(Float64, size(array.positions, 2), 2)
    for (n, p) in enumerate(eachcol(array.positions))
        J[n, 1] = 2*π * dot(p, du)
        J[n, 2] = 2*π * dot(p, dv)
    end
    return J
end

steering_derivative_norm(array, point::UVDirection) = opnorm(steering_derivative_matrix(array, point))

function pointing_delta(tol::Tolerances, array, point::ThetaDirection)
    return tol.pointing[1] * steering_derivative_norm(physical_array(array), point)
end

function pointing_delta(tol::Tolerances, array, point::UVDirection)
    J = steering_derivative_matrix(physical_array(array), point)
    return opnorm(J * Diagonal(collect(tol.pointing)))
end

function channel_delta(tol::Tolerances, N::Integer)
    c = hypot(tol.amplitude, tol.phase)
    return c * sqrt(N)
end

function position_delta(tol::Tolerances, point, N::Integer)
    u, v, w = direction_vector(point)
    ρx, ρy, ρz = tol.position
    return 2*π * sqrt(ρx^2 * u^2 + ρy^2 * v^2 + ρz^2 * w^2) * sqrt(N)
end

function combine_deltas(tol::Tolerances, deltas)
    tol.combine === :sum && return sum(deltas)
    return sqrt(sum(abs2, deltas))
end

function delta_at(tol::Tolerances, array, point)
    physical = physical_array(array)
    N = size(physical.positions, 2)
    return combine_deltas(tol, (
        pointing_delta(tol, physical, point),
        channel_delta(tol, N),
        position_delta(tol, point, N),
    ))
end

robust_norm_entries(array, vars::SparseVariables) = robust_norm_entries(array, vars.variables)

function robust_norm_entries(array::ArrayGeometry, vars::WeightVariables)
    entries = AffExpr[]
    append!(entries, 1.0 .* vars.w_re)
    append!(entries, 1.0 .* vars.w_im)
    return entries
end

function robust_norm_entries(array::ArrayGeometry, vars::AmplitudeVariables)
    return 1.0 .* vars.a
end

function robust_norm_entries(array::SymmetricArray, vars::WeightVariables)
    entries = AffExpr[]
    for n in eachindex(vars.w_re)
        scale = is_origin(view(array.positions, :, n)) ? 1.0 : sqrt(2.0)
        push!(entries, scale * vars.w_re[n])
    end
    for n in eachindex(vars.w_im)
        scale = is_origin(view(array.positions, :, n)) ? 1.0 : sqrt(2.0)
        push!(entries, scale * vars.w_im[n])
    end
    return entries
end

function robust_norm_entries(array::SymmetricArray, vars::AmplitudeVariables)
    entries = AffExpr[]
    for n in eachindex(vars.a)
        scale = is_origin(view(array.positions, :, n)) ? 1.0 : sqrt(2.0)
        push!(entries, scale * vars.a[n])
    end
    return entries
end

robust_margin!(model, array, weights, vars, formulation, ::Nothing) = nothing

function robust_margin!(model, array, weights, vars, formulation, tol::Tolerances)
    formulation isa SOCP || error("robustness requires SOCP formulation.")
    τ = @variable(model, lower_bound = 0.0)
    entries = robust_norm_entries(array, vars)
    @constraint(model, [1.0 * τ; entries] in SecondOrderCone())
    return RobustMargin(tol, τ)
end

robust_bound(bound, ::Nothing, array, point) = bound

function robust_bound(bound, margin::RobustMargin, array, point)
    δ = delta_at(margin.tolerances, array, point)
    return bound - δ * margin.norm_bound
end
