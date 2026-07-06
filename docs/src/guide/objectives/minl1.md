# MinL1

Minimize the L1 norm of the optimization variables.

## Syntax

```julia
MinL1(; sum_limit = nothing)
```

## Arguments

| Argument | Description |
|---|---|
| `sum_limit` | Optional upper bound on the sum of the L1 auxiliary variables. |

## Description

Introduces one auxiliary variable per optimization variable and minimizes their sum.
This promotes sparse or low-amplitude excitations while keeping the problem linear
when used with [`LP`](@ref "LP").

`MinL1` is a direct objective. For stronger sparsity promotion through repeated
reweighting, use [`IterativeReweightedL1`](@ref "IterativeReweightedL1").

## Example

```julia
using ArraySynthesis; using ArraySynthesis: °, dB; using HiGHS

array = uniform_linear_array(32, d = 0.5)
p     = pattern(beam(0°),
                sidelobes(region(-90°..(-5°), 1°), -25dB),
                sidelobes(region(5°..90°,     1°), -25dB))

result = synthesize(array, p, MinL1(), ComplexWeights(), LP(), HiGHS.Optimizer)
```

## Related

**Sparse objectives:** [`MinWeightedL1`](@ref "MinWeightedL1"),
[`IterativeReweightedL1`](@ref "IterativeReweightedL1")

**Compatible formulation:** [`LP`](@ref "LP")
