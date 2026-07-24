# MinL1

Minimize a convex L1-type bound on the optimization variables.

## Syntax

```julia
MinL1(sum_limit)
MinL1(; sum_limit = nothing)
```

## Arguments

| Argument | Description |
|---|---|
| `sum_limit` | Optional upper bound on the sum of the L1 auxiliary variables. |

## Description

Introduces one auxiliary variable per optimization variable and minimizes their sum. With [`LP`](@ref "LP") and [`QP`](@ref "QP"), these auxiliaries are imposed with linear inequalities. With [`SOCP`](@ref "SOCP"), complex pairs can be bounded with a second-order cone.

When `sum_limit` is provided, the model also imposes an upper bound on `sum(t)`.

`MinL1` is a direct convex objective. It is not the sparse/thinned-array strategy. For active-element thinning through repeated reweighting, use [`IterativeReweightedL1`](@ref "IterativeReweightedL1").

## Example

```julia
using ArraySynthesis; using ArraySynthesis: °, dB; using HiGHS

array = uniform_linear_array(32, d = 0.5)
p = pattern(beam(0°),
            sidelobes(region(-90°..(-5°), 1°), -25dB),
            sidelobes(region(5°..90°,     1°), -25dB))

result = synthesize(array, p, MinL1(), ComplexWeights(), LP(), HiGHS.Optimizer)
```