struct SynthesisResult
    weights::Vector{ComplexF64}
    status::MOI.TerminationStatusCode
    objective_value::Float64
    model::Model
end

struct IterativeSynthesisResult
    weights::Vector{ComplexF64}
    converged::Bool
    iterations::Int
end


function extract_phase_direction(pattern::Pattern)
    # Try main beam direction
    if !isempty(pattern.beams)
        return pattern.beams[1].direction
    end
    # Try shaped beam region center (average of all points in region)
    if !isempty(pattern.shaped_beams)
        pts = pattern.shaped_beams[1].region.points
        if pts[1] isa ThetaDirection
            avg_θ = sum(p.θ for p in pts) / length(pts)
            return ThetaDirection(avg_θ)
        elseif pts[1] isa UVDirection
            avg_u = sum(p.u for p in pts) / length(pts)
            avg_v = sum(p.v for p in pts) / length(pts)
            return UVDirection(avg_u, avg_v)
        end
    end
    return ThetaDirection(0.0)
end

function resolve_phase_reference(weights::ProgressivePhaseAmplitude, pattern::Pattern)
    weights.β === nothing && return ProgressivePhaseAmplitude(extract_phase_direction(pattern))
    return weights
end

resolve_phase_reference(weights, pattern) = weights

struct MultiPatternResult
    weights::Vector{Vector{ComplexF64}}
    converged::Bool
    iterations::Int
end


function extract_weights(vars::WeightVariables)
    re = value.(vars.w_re)
    isempty(vars.w_im) && return ComplexF64.(re)
    return complex.(re, value.(vars.w_im))
end

extract_weights(vars::AmplitudeVariables) = ComplexF64.(value.(vars.a))
extract_weights(vars::SparseVariables) = extract_weights(vars.variables)


function filter_beam_objectives(pattern, objective)
    skip = Set(objective_beam_directions(objective, pattern))
    isempty(skip) && return pattern
    beams = filter(b -> b.direction ∉ skip, pattern.beams)
    return Pattern(beams, pattern.shaped_beams, pattern.null_directions, pattern.sidelobe_regions)
end

function solve_model(array, pattern, objective::DirectObjective, weights, formulation, solver; solver_options = nothing)
    check_formulation(objective, formulation)
    model = Model(solver)
    if !isnothing(solver_options)
        for (key, value) in solver_options
            set_optimizer_attribute(model, string(key), value)
        end
    end
    set_silent(model)
    vars = variables!(model, array, weights, formulation)
    constraints!(model, filter_beam_objectives(pattern, objective), array, weights, vars, formulation)
    aux = objective!(model, objective, pattern, array, weights, vars, formulation)
    optimize!(model)
    return model, vars, aux
end


function synthesize(array, pattern, objective::DirectObjective, weights, formulation, solver; solver_options = nothing)
    weights = resolve_phase_reference(weights, pattern)
    model, vars, _ = solve_model(array, pattern, objective, weights, formulation, solver; solver_options = solver_options)
    is_solved_and_feasible(model) || @warn "Synthesis failed: $(termination_status(model))"
    return SynthesisResult(extract_weights(vars), termination_status(model), objective_value(model), model)
end

function synthesize(array, pattern, objective::IterativeReweightedL1, weights, formulation, solver; solver_options = nothing)
    weights = resolve_phase_reference(weights, pattern)
    N = size(array.positions, 2)
    alpha = ones(N)
    delta = 0.0
    n_active_prev = -1
    n_stable = 0
    w = zeros(ComplexF64, N)

    for it in 1:objective.max_iter
        model, vars, t_vars = solve_model(array, pattern, MinWeightedL1(alpha), weights, formulation, solver; solver_options = solver_options)
        is_solved_and_feasible(model) || @warn "IterativeReweightedL1 failed at iteration $it: $(termination_status(model))"

        t = value.(t_vars)
        w = extract_weights(vars)
        it == 1 && (delta = objective.d0 * maximum(t))
        alpha = 1.0 ./ (t .+ delta)

        n_active = count(t .>= maximum(t) / 1e5)
        n_stable = n_active == n_active_prev ? n_stable + 1 : 0
        n_active_prev = n_active
        n_stable >= objective.stable_iterations && return IterativeSynthesisResult(w, true, it)
    end

    return IterativeSynthesisResult(w, false, objective.max_iter)
end


function array_factor(array::ArrayGeometry, ::Union{ComplexWeights, RealWeights}, w, pts)
    return steering_matrix(array, direction_matrix(pts)) * w
end

function array_factor(array::ArrayGeometry, weights::ProgressivePhaseAmplitude, w, pts)
    A = steering_matrix(array, direction_matrix(pts) .- phase_direction(weights.β))
    return A * real.(w)
end

function array_factor(array::SymmetricArray, ::ConjugateSymmetricWeights, w, pts)
    A_cos, A_sin = steering_matrix(array, direction_matrix(pts))
    return [sum(A_cos[m, n] * real(w[n]) + A_sin[m, n] * imag(w[n]) for n in axes(A_cos, 2)) for m in eachindex(pts)]
end

function array_factor(array::SymmetricArray, weights::ProgressivePhaseAmplitude, w, pts)
    A_cos, _ = steering_matrix(array, direction_matrix(pts) .- phase_direction(weights.β))
    a = real.(w)
    return [sum(A_cos[m, n] * a[n] for n in axes(A_cos, 2)) for m in eachindex(pts)]
end

# synthesize: IterativePatternLeastSquares
# Each iteration: solve MinFieldError with complex target, then update target phase
# from the current AF (Gerchberg-Saxton-style phase update).
function synthesize(array, pattern, objective::IterativePatternLeastSquares, weights, formulation, solver; solver_options = nothing)
    weights = resolve_phase_reference(weights, pattern)
    pts = objective.region.points
    mag = abs.(objective.target)
    complex_target = ComplexF64.(objective.target)
    w = zeros(ComplexF64, size(array.positions, 2))

    for it in 1:objective.max_iter
        model, vars, _ = solve_model(array, pattern, MinFieldError(pts, complex_target), weights, formulation, solver; solver_options = solver_options)
        is_solved_and_feasible(model) ||
            error("IterativePatternLeastSquares failed at iteration $it: $(termination_status(model))")

        w_new = extract_weights(vars)
        af = array_factor(array, weights, w_new, pts)
        complex_target = mag .* cis.(angle.(af))

        norm(w_new - w) < 1e-6 * (1 + norm(w_new)) && return IterativeSynthesisResult(w_new, true, it)
        w = w_new
    end

    return IterativeSynthesisResult(w, false, objective.max_iter)
end

# IterativeFloorSynthesis helpers
# Orchard-Elliott-Stern peak correction. Only supported for ArrayGeometry (full array).
function slope_row_theta(array::ArrayGeometry, θ, A0_row)
    N = size(array.positions, 2)
    return [-im * 2π * (cos(θ) * array.positions[1, n] - sin(θ) * array.positions[3, n]) * A0_row[n] for n in 1:N]
end

function local_sidelobe_peaks(mag, M; left_endpoint = true, right_endpoint = true, score = mag)
    peaks = Int[]
    length(mag) >= 2 && left_endpoint && mag[1] > mag[2] && push!(peaks, 1)
    for i in 2:length(mag)-1
        if mag[i] > mag[i-1] && mag[i] > mag[i+1]
            push!(peaks, i)
        end
    end
    length(mag) >= 2 && right_endpoint && mag[end] > mag[end-1] && push!(peaks, length(mag))
    isempty(peaks) && return Int[]
    order = sortperm(score[peaks], rev=true)
    return peaks[order[1:min(M, length(order))]]
end

function sidelobe_peak_indices(regions, levels, af, M)
    peaks = Int[]
    offset = 0
    for sl in regions
        n = length(sl.region.points)
        inds = (offset + 1):(offset + n)
        pts = sl.region.points
        left_endpoint = first(pts) isa ThetaDirection && first(pts).θ <= -π/2 + 1e-12
        right_endpoint = last(pts) isa ThetaDirection && last(pts).θ >= π/2 - 1e-12
        local_peaks = local_sidelobe_peaks(abs.(af[inds]), n;
            left_endpoint = left_endpoint,
            right_endpoint = right_endpoint,
            score = abs.(af[inds]) ./ levels[inds])
        append!(peaks, offset .+ local_peaks)
        offset += n
    end
    isempty(peaks) && return Int[]
    order = sortperm(abs.(af[peaks]) ./ levels[peaks], rev=true)
    return peaks[order[1:min(M, length(order))]]
end

function solve_floor_energy_qp(A_side, A0_row, dA0_row, null_rows, solver)
    N = size(A_side, 2)
    model = Model(solver)
    set_silent(model)
    w_re = @variable(model, [1:N])
    w_im = @variable(model, [1:N])

    af_re(row) = sum(real(row[n]) * w_re[n] - imag(row[n]) * w_im[n] for n in 1:N)
    af_im(row) = sum(real(row[n]) * w_im[n] + imag(row[n]) * w_re[n] for n in 1:N)

    @constraint(model, af_re(A0_row) == 1.0)
    @constraint(model, af_im(A0_row) == 0.0)
    @constraint(model, af_re(dA0_row) == 0.0)

    for row in null_rows
        @constraint(model, af_re(row) == 0.0)
        @constraint(model, af_im(row) == 0.0)
    end

    @objective(model, Min, sum(af_re(A_side[m, :])^2 + af_im(A_side[m, :])^2 for m in axes(A_side, 1)))

    optimize!(model)
    is_solved_and_feasible(model) || @warn "IterativeFloorSynthesis initialization failed: $(termination_status(model))"
    return value.(w_re) .+ im .* value.(w_im)
end

function solve_floor_correction_qp(A_side, peak_rows, peak_targets, A0_row, dA0_row, null_rows, solver)
    N = size(A_side, 2)
    model = Model(solver)
    set_silent(model)
    dw_re = @variable(model, [1:N])
    dw_im = @variable(model, [1:N])

    row_re(row) = sum(real(row[n]) * dw_re[n] - imag(row[n]) * dw_im[n] for n in 1:N)
    row_im(row) = sum(real(row[n]) * dw_im[n] + imag(row[n]) * dw_re[n] for n in 1:N)

    for k in axes(peak_rows, 1)
        row = peak_rows[k, :]
        target = peak_targets[k]
        @constraint(model, row_re(row) == real(target))
        @constraint(model, row_im(row) == imag(target))
    end

    @constraint(model, row_re(A0_row) == 0.0)
    @constraint(model, row_im(A0_row) == 0.0)
    @constraint(model, row_re(dA0_row) == 0.0)

    for row in null_rows
        @constraint(model, row_re(row) == 0.0)
        @constraint(model, row_im(row) == 0.0)
    end

    residual_re(m) = row_re(A_side[m, :])
    residual_im(m) = row_im(A_side[m, :])
    @objective(model, Min, sum(residual_re(m)^2 + residual_im(m)^2 for m in axes(A_side, 1)))

    optimize!(model)
    is_solved_and_feasible(model) || return nothing
    return value.(dw_re) .+ im .* value.(dw_im)
end

function synthesize(array::SymmetricArray, pattern, objective::IterativeFloorSynthesis, weights, formulation, solver)
    error("IterativeFloorSynthesis requires ArrayGeometry. Use materialize(array) to get the full array.")
end

function synthesize(array::ArrayGeometry, pattern, objective::IterativeFloorSynthesis, weights, formulation, solver)
    weights = resolve_phase_reference(weights, pattern)
    isempty(pattern.sidelobe_regions) && error("IterativeFloorSynthesis requires at least one sidelobe region in the pattern.")
    objective.direction isa ThetaDirection || error("IterativeFloorSynthesis only supports ThetaDirection.")

    side_pts = reduce(vcat, [sl.region.points for sl in pattern.sidelobe_regions])
    side_levels = reduce(vcat, [sl.upper for sl in pattern.sidelobe_regions])
    null_pts = [n.direction for n in pattern.null_directions]

    N = size(array.positions, 2)
    U_side = direction_matrix(side_pts)
    A_side = steering_matrix(array, U_side)
    A0_row = vec(steering_matrix(array, direction_matrix([objective.direction])))
    dA0_row = slope_row_theta(array, objective.direction.θ, A0_row)
    null_rows = [vec(steering_matrix(array, direction_matrix([p]))) for p in null_pts]

    w = solve_floor_energy_qp(A_side, A0_row, dA0_row, null_rows, solver)
    M = max(N - 2 - length(null_pts), 1)

    for it in 1:objective.max_iter
        af = A_side * w
        maximum(abs.(af) ./ side_levels) <= 1 + objective.tol && return IterativeSynthesisResult(w, true, it)

        peaks = sidelobe_peak_indices(pattern.sidelobe_regions, side_levels, af, M)
        isempty(peaks) && return IterativeSynthesisResult(w, true, it)

        dW = nothing
        for npeaks in length(peaks):-1:1
            selected = peaks[1:npeaks]
            ci = af[selected]
            target = side_levels[selected]
            fi = (target .- abs.(ci)) .* ci ./ abs.(ci)
            dW = solve_floor_correction_qp( A_side, A_side[selected, :], fi, A0_row, dA0_row, null_rows, solver)
            dW === nothing || break
        end
        dW === nothing && return IterativeSynthesisResult(w, false, it)
        w .+= dW

        norm(dW) < objective.tol && return IterativeSynthesisResult(w, true, it)
    end

    return IterativeSynthesisResult(w, false, objective.max_iter)
end

# synthesize: MultiPatternReweightedL1
# K patterns share one L1 auxiliary variable t. Each pattern has independent
# weights, but an element is counted through the largest excitation it needs.
function solve_multipattern_model(array, patterns, weights, formulation, solver, alpha)
    model = Model(solver)
    set_silent(model)
    N = size(array.positions, 2)
    t = @variable(model, [1:N])

    all_vars = map(patterns) do pat
        vars = variables!(model, array, weights, formulation)
        constraints!(model, pat, array, weights, vars, formulation)
        l1_bound!(model, t, vars, formulation)
        vars
    end

    @objective(model, Min, sum(alpha[n] * t[n] for n in eachindex(t)))
    optimize!(model)
    return model, all_vars, value.(t)
end

function synthesize(array, patterns::AbstractVector{<:Pattern}, objective::MultiPatternReweightedL1, weights, formulation, solver)
    # Resolve phase reference from first pattern
    weights = resolve_phase_reference(weights, patterns[1])
    N = size(array.positions, 2)
    alpha = ones(N)
    delta = 0.0
    n_active_prev = -1
    n_stable = 0
    all_weights = [zeros(ComplexF64, N) for _ in patterns]

    for it in 1:objective.max_iter
        model, all_vars, t = solve_multipattern_model(array, patterns, weights, formulation, solver, alpha)
        is_solved_and_feasible(model) ||
            error("MultiPatternReweightedL1 failed at iteration $it: $(termination_status(model))")

        all_weights = [extract_weights(v) for v in all_vars]
        it == 1 && (delta = objective.d0 * maximum(t))
        alpha = 1.0 ./ (t .+ delta)

        n_active = count(t .>= maximum(t) / 1e5)
        n_stable = n_active == n_active_prev ? n_stable + 1 : 0
        n_active_prev = n_active
        n_stable >= objective.stable_iterations && return MultiPatternResult(all_weights, true, it)
    end

    return MultiPatternResult(all_weights, false, objective.max_iter)
end
