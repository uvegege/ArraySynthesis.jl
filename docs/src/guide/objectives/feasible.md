# Feasible

Find any feasible solution without optimizing any objective.

## Syntax

```julia
Feasible()
```

## Description

Sets the objective to zero and returns the first feasible point found by the solver.
Useful to check whether the pattern constraints can be satisfied at all, or to obtain
a baseline excitation before switching to a more specific objective.

## Example

```julia
using ArraySynthesis; using ArraySynthesis: °, dB; using HiGHS

array = uniform_linear_array(32, d = 0.5)
p     = pattern(beam(0°),
                sidelobes(region(-90°..(-5°), 1°), -25dB),
                sidelobes(region(5°..90°,     1°), -25dB))

result = synthesize(array, p, Feasible(), ComplexWeights(), LP(), HiGHS.Optimizer)
```

## Related

**Alternative objectives:** [`MaxGain`](@ref "MaxGain"), [`MinSLL`](@ref "MinSLL"),
[`MinPower`](@ref "MinPower")
