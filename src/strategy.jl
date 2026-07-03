abstract type SynthesisMethod <: AbstractObjective end

objective!(model, objective::SynthesisMethod, pattern, array, weights, vars, formulation) = error("$(typeof(objective)) requires a synthesis loop, not a single compiled model.")

# Reweighted L1 minimization for sparse arrays (Candès-Wakin-Boyd).
# Iteratively reweights excitation activities to push small elements toward zero.
struct IterativeReweightedL1{D} <: SynthesisMethod
    max_iter::Int
    d0::D
    stable_iterations::Int
end

function IterativeReweightedL1(; max_iter = 15, d0 = 1e-3, stable_iterations = 3)
   return IterativeReweightedL1(max_iter, d0, stable_iterations)
end

# Iterative LP/QP for shaped beam synthesis (Fuchs method).
# Each iteration updates the phase reference from the previous solution.
struct IterativePatternLeastSquares{R, T} <: SynthesisMethod
    region::R
    target::Vector{T}
    max_iter::Int
end

function IterativePatternLeastSquares(region::Region, target::AbstractVector; max_iter::Int = 15)
    length(target) == length(region.points) || error("Target length must match region points.")
    return IterativePatternLeastSquares(region, collect(target), max_iter)
end

function IterativePatternLeastSquares(region::Region, target::Number; max_iter = 15) 
    return IterativePatternLeastSquares(region, fill(target, length(region.points)); max_iter)
end

# Orchard-Elliott-Stern iterative peak correction.
# Imposes a real SLL floor without a phase reference by correcting the highest
# sidelobe peaks one QP perturbation at a time.
struct IterativeFloorSynthesis{P} <: SynthesisMethod
    direction::P
    sll::Float64
    max_iter::Int
    tol::Float64
    stable_iterations::Int
end

function IterativeFloorSynthesis(dir; sll = -30.0dB, max_iter::Int = 15, tol::Real = 1e-4, stable_iterations::Int = 3)
    return IterativeFloorSynthesis(direction(dir), Float64(sll), max_iter, Float64(tol), stable_iterations)
end

# Reweighted L1 for K patterns sharing the same active element support.
# One shared L1 auxiliary variable, independent weight variables per pattern.
struct MultiPatternReweightedL1{D} <: SynthesisMethod
    max_iter::Int
    d0::D
    stable_iterations::Int
end

function MultiPatternReweightedL1(; max_iter::Int = 15, d0 = 1e-3, stable_iterations = 3)
    return MultiPatternReweightedL1(max_iter, d0, stable_iterations)
end
