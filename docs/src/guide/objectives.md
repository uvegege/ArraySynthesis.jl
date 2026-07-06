# Objectives

The `AbstractObjective` determines what the solver minimizes or maximizes, subject to
the constraints defined by [`Pattern`](@ref). Pass it as the third argument to
[`synthesize`](@ref).

!!! note "Provisional guide"
    This page is a selection guide. The individual objective pages contain the
    syntax, arguments, examples, and implementation notes.

## Choosing an Objective

| Goal | Objective | Typical formulation |
|---|---|---|
| Check feasibility only | [`Feasible`](@ref "Feasible") | `LP` |
| Minimize peak sidelobe level | [`MinSLL`](@ref "MinSLL") | `LP` |
| Maximize field at a look direction | [`MaxGain`](@ref "MaxGain") | `LP` |
| Minimize the L1 norm of the weights | [`MinL1`](@ref "MinL1") | `LP` |
| Minimize a weighted L1 norm | [`MinWeightedL1`](@ref "MinWeightedL1") | `LP` |
| Sparse array synthesis | [`IterativeReweightedL1`](@ref "IterativeReweightedL1") | `LP` |
| Sparse synthesis for multiple patterns | [`MultiPatternReweightedL1`](@ref "MultiPatternReweightedL1") | `LP` |
| Minimize excitation power | [`MinPower`](@ref "MinPower") | `QP` |
| Minimize radiated power over a region | [`MinIntegratedPower`](@ref "MinIntegratedPower") | `QP` |
| Fit a complex field reference | [`MinFieldError`](@ref "MinFieldError") | `QP` |
| Iterative shaped-beam least squares | [`IterativePatternLeastSquares`](@ref "IterativePatternLeastSquares") | `LP` |
| Iterative low-sidelobe synthesis | [`IterativeFloorSynthesis`](@ref "IterativeFloorSynthesis") | `QP` |
| Maximize directivity | [`MaxDirectivity`](@ref "MaxDirectivity") | `SOCP` |

## Direct Objectives

Direct objectives solve one convex optimization problem.

Use [`Feasible`](@ref "Feasible") when the pattern constraints are the result and no
secondary criterion is needed. Use [`MinSLL`](@ref "MinSLL") when sidelobe level is the
quantity to optimize rather than a fixed constraint. Use [`MaxGain`](@ref "MaxGain")
when the main concern is gain at one or more look directions.

Use [`MinL1`](@ref "MinL1") and [`MinWeightedL1`](@ref "MinWeightedL1") when the direct
objective is sparse or low-amplitude weights. The iterative sparse objectives build on
the same idea when a single L1 solve is not selective enough.

Quadratic objectives such as [`MinPower`](@ref "MinPower"),
[`MinIntegratedPower`](@ref "MinIntegratedPower"), and [`MinFieldError`](@ref "MinFieldError")
need at least [`QP`](@ref "QP"). [`MaxDirectivity`](@ref "MaxDirectivity") is formulated
as a conic problem and requires [`SOCP`](@ref "SOCP").

## Iterative Objectives

Iterative objectives solve a sequence of convex problems. They are useful when the
desired behavior is not represented by a single convex objective, such as shaped-beam
phase recovery or sparse element selection.

[`IterativePatternLeastSquares`](@ref "IterativePatternLeastSquares") updates the phase
reference of a shaped beam. [`IterativeReweightedL1`](@ref "IterativeReweightedL1") and
[`MultiPatternReweightedL1`](@ref "MultiPatternReweightedL1") promote sparse active
element sets. [`IterativeFloorSynthesis`](@ref "IterativeFloorSynthesis") follows a
floor-constraint strategy for narrow-beam low-sidelobe synthesis.

## Related

**Problem setup:** [`Pattern`](@ref), [`synthesize`](@ref)

**Model choices:** [`Excitations`](@ref "Excitations"), [`Formulations`](@ref "Formulations")
