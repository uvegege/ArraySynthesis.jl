# MinWeightedL1

Minimize a weighted convex L1-type bound on the optimization variables.

## Syntax

```julia
MinWeightedL1(alpha; sum_limit = nothing)
```

## Arguments

| Argument | Description |
|---|---|
| `alpha` | Weight for each optimization variable. Its length must match the number of variables. |
| `sum_limit` | Optional upper bound on the sum of the L1 auxiliary variables. |

## Description

`MinWeightedL1` is the weighted version of [`MinL1`](@ref "MinL1"). It introduces the same auxiliary bounds and minimizes `sum(alpha[n] * t[n])`, so elements or variables can be penalized more heavily than others.

When `sum_limit` is provided, the model also imposes an upper bound on `sum(t)`.

This is a direct convex objective, not the sparse/thinned-array strategy. Sparse thinning is implemented by [`IterativeReweightedL1`](@ref "IterativeReweightedL1") and [`MultiPatternReweightedL1`](@ref "MultiPatternReweightedL1"), which repeatedly update these weights.

## Example

```julia
using ArraySynthesis; using ArraySynthesis: °, dB; using HiGHS

array = uniform_linear_array(16, d = 0.5)
p = pattern(beam(0°),
            sidelobes(region(-90°..(-5°), 1°), -25dB),
            sidelobes(region(5°..90°, 1°), -25dB))
alpha = ones(2 * 16)

result = synthesize(array, p, MinWeightedL1(alpha), ComplexWeights(), LP(), HiGHS.Optimizer)
```
