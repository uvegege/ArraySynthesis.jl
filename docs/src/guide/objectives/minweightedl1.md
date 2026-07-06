# MinWeightedL1

Minimize a weighted L1 norm of the optimization variables.

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

`MinWeightedL1` is the weighted version of [`MinL1`](@ref "MinL1"). It is useful when
some elements or variables should be penalized more heavily than others.

The objective remains linear when used with [`LP`](@ref "LP").

## Example

```julia
using ArraySynthesis; using ArraySynthesis: °, dB; using HiGHS

array = uniform_linear_array(16, d = 0.5)
p     = pattern(beam(0°),
                sidelobes(region(-90°..(-5°), 1°), -25dB),
                sidelobes(region(5°..90°,     1°), -25dB))
alpha = ones(2 * 16)

result = synthesize(array, p, MinWeightedL1(alpha), ComplexWeights(), LP(), HiGHS.Optimizer)
```

## Related

**Sparse objectives:** [`MinL1`](@ref "MinL1"),
[`IterativeReweightedL1`](@ref "IterativeReweightedL1")

**Compatible formulation:** [`LP`](@ref "LP")
