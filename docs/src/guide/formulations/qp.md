# QP

Quadratic programming formulation with polyhedral modulus constraint.

## Syntax

```julia
QP()
QP(polygon_faces = F)
```

## Arguments

| Argument | Description |
|---|---|
| `polygon_faces` | Number of polygon faces for the modulus approximation (default `8`). |

## Description

Uses the same polyhedral (linear) modulus constraint as `LP`. The difference is the
objective: `QP` allows a quadratic term, enabling objectives such as `MinPower`,
`MinIntegratedPower`, and `MinFieldError`. The constraints themselves are always linear.

Required by `IterativeFloorSynthesis`.

Suggested solvers: `HiGHS.Optimizer`, `Mosek.Optimizer`

## Example

```julia
using ArraySynthesis; using ArraySynthesis: °, dB; using HiGHS

array = uniform_linear_array(32, d = 0.5)
p     = pattern(beam(0°),
                sidelobes(region(-90°..(-5°), 1°), -20dB),
                sidelobes(region(5°..90°,     1°), -20dB),
                null(-40°, level = -55dB))

result = synthesize(array, p, MinPower(), ComplexWeights(), QP(), HiGHS.Optimizer)
```

## Related

**Alternative formulations:** [`LP`](@ref "LP"), [`SOCP`](@ref "SOCP")

**Objective requiring QP:** [`IterativeFloorSynthesis`](@ref "IterativeFloorSynthesis")

**Background:** [Modulus Constraints per Formulation](@ref)
