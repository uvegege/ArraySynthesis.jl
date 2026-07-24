# Objectives

The `AbstractObjective` determines what the solver minimizes or maximizes, subject to the constraints defined by [`Pattern`](@ref). Pass it as the third argument to [`synthesize`](@ref).

Use this page to choose the objective. For exact constructor signatures and argument details, see the [API Reference](@ref).

## Choosing an Objective

| Goal | Objective | Typical formulation |
|---|---|---|
| Check feasibility only | [`Feasible`](@ref "Feasible") | `LP` |
| Minimize peak sidelobe level | [`MinSLL`](@ref "MinSLL") | `LP` |
| Maximize array-factor response at a look direction | [`MaxAF`](@ref "MaxAF") | `LP` |
| Minimize an L1-type excitation bound | [`MinL1`](@ref "MinL1") | `LP` |
| Minimize a weighted L1-type excitation bound | [`MinWeightedL1`](@ref "MinWeightedL1") | `LP` |
| Sparse array synthesis | [`IterativeReweightedL1`](@ref "IterativeReweightedL1") | `LP` |
| Sparse synthesis for multiple patterns | [`MultiPatternReweightedL1`](@ref "MultiPatternReweightedL1") | `LP` |
| Minimize excitation power | [`MinPower`](@ref "MinPower") | `QP` |
| Minimize array-factor power over a region | [`MinIntegratedPower`](@ref "MinIntegratedPower") | `QP` |
| Fit a complex array-factor reference | [`MinFieldError`](@ref "MinFieldError") | `QP` |
| Iterative low-sidelobe synthesis | [`IterativeFloorSynthesis`](@ref "IterativeFloorSynthesis") | `QP` |

## Direct Objectives

Direct objectives solve one convex optimization problem.

Use [`Feasible`](@ref "Feasible") when the pattern constraints are the result and no secondary criterion is needed. Use [`MinSLL`](@ref "MinSLL") when sidelobe level is the quantity to optimize rather than a fixed constraint. Use [`MaxAF`](@ref "MaxAF") when the main concern is gain at one or more look directions.

Use [`MinL1`](@ref "MinL1") and [`MinWeightedL1`](@ref "MinWeightedL1") when the direct objective is a convex bound on excitation variables. They do not implement thinning or active-element selection by themselves; sparse/thinned synthesis is handled by the iterative reweighted strategies.

Quadratic objectives such as [`MinPower`](@ref "MinPower"), [`MinIntegratedPower`](@ref "MinIntegratedPower"), and [`MinFieldError`](@ref "MinFieldError") need at least [`QP`](@ref "QP").

### Constraint satisfaction

[`Feasible`](@ref "Feasible") only asks whether the requested pattern constraints can be satisfied. It is useful when another criterion is unnecessary or when debugging a new mask.

### Gain and sidelobe level

[`MaxAF`](@ref "MaxAF") maximizes the real array-factor response in one or more look directions. [`MinSLL`](@ref "MinSLL") introduces auxiliary sidelobe-level variables over one or more regions and minimizes their sum.

### Norm and power objectives

[`MinL1`](@ref "MinL1") and [`MinWeightedL1`](@ref "MinWeightedL1") minimize sums of
auxiliary variables that bound excitation magnitudes or components. [`MinPower`](@ref "MinPower") minimizes excitation power, while [`MinIntegratedPower`](@ref "MinIntegratedPower") minimizes array-factor power over sampled directions.

### Array-factor fitting

[`MinFieldError`](@ref "MinFieldError") fits a scalar, vector, or complex array-factor reference over a region.

## Iterative Objectives

Iterative objectives solve a sequence of convex problems. They are useful when the desired behavior is not represented by a single convex objective, such as sparse element selection or peak-by-peak sidelobe correction.

[`IterativeReweightedL1`](@ref "IterativeReweightedL1") and [`MultiPatternReweightedL1`](@ref "MultiPatternReweightedL1") promote sparse active element sets. [`IterativeFloorSynthesis`](@ref "IterativeFloorSynthesis") follows a floor-constraint strategy for narrow-beam low-sidelobe synthesis.
