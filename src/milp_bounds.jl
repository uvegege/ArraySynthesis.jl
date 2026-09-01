function collect_bounded_pattern(array, pattern::Pattern{P,T}, objective = nothing) where {P,T}
    points = P[]
    bounds = T[]

    for beam in pattern.beams
        push!(points, beam.direction)
        push!(bounds, beam.gain)
    end

    for beam in pattern.shaped_beams
        append!(points, beam.region.points)
        append!(bounds, map(beam.target) do target
            b1 = target * beam.ripple
            b2 = target / beam.ripple
            max(abs(b1), abs(b2))
        end)
    end

    for sidelobes in pattern.sidelobe_regions
        append!(points, sidelobes.region.points)
        append!(bounds, sidelobes.upper)
    end

    if objective isa MinSLL && objective.upper_bound !== nothing
        for region in objective.regions
            append!(points, region.points)
            append!(bounds, fill(objective.upper_bound, length(region.points)))
        end
    end

    directions = direction_matrix(points)
    A_cos, A_sin = steering_matrix(array, directions)

    N_rep = size(array.positions, 2)
    A_full = hcat(A_cos, A_sin)

    return A_full, bounds, N_rep, length(points)
end

function analytic_bounds(A, b, N_rep; rtol = sqrt(eps(Float64)))
    m, n = size(A)

    n == 2 * N_rep || throw(DimensionMismatch("Expected $((2 * N_rep)) columns, got $n."))
    length(b) == m || throw(DimensionMismatch("Expected $m bounds, got $(length(b))."))
    isempty(b) && return nothing

    s = svdvals(A)
    length(s) == n || return nothing

    σmax = first(s)
    σmin = last(s)
    σmin > rtol * σmax || return nothing

    U_norm = norm(b) / σmin # ||x||₂ <= ||b||₂ / σmin(A)

    # x = A† y, |y| <= b
    # |x| <= |A†| b
    x_bound = abs.(pinv(A; rtol)) * abs.(b)

    re_bound = @view x_bound[1:N_rep]
    im_bound = @view x_bound[(N_rep + 1):(2*N_rep)]

    U = hypot.(re_bound, im_bound)

    return (per_element = U, global_bound = maximum(U), 
    norm_bound = U_norm, sigma_min = σmin, 
    sigma_max = σmax, condition = σmax / σmin)
end

function amplitude_upper_bounds(array, pattern, objective = nothing; kwargs...)
    A, b, Nrep, _ = collect_bounded_pattern(array, pattern, objective)
    return analytic_bounds(A, b, Nrep; kwargs...)
end