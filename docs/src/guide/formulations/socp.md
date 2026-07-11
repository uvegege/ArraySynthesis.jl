# SOCP

Second-order cone programming formulation with exact modulus constraint.

## Syntax

```julia
SOCP()
```

## Description

Imposes $|AF(\hat{r}_p)| \leq t_p$ exactly as a second-order cone constraint:

```math
\|(AF_{re}(\hat{r}_p),\, AF_{im}(\hat{r}_p))\|_2 \leq t_p
```

No polygon approximation; the modulus bound is always exact.

Required by `MaxDirectivity`. Compatible with linear and quadratic objectives.
Typically slower than `LP`/`QP` on large models but more accurate for complex AF.

Also required when using [`robust`](@ref "Robustness"), because robust bounds use
one shared norm cone for the excitation vector.

Suggested solvers: `Clarabel.Optimizer`, `Mosek.Optimizer`

## Example

```julia
using ArraySynthesis; using ArraySynthesis: °, dB; using Clarabel

array = uniform_linear_array(32, d = 0.5)
p     = pattern(beam(0°),
                sidelobes(region(-90°..(-5°), 1°), -20dB),
                sidelobes(region(5°..90°,     1°), -20dB),
                null(-40°, level = -55dB))

result = synthesize(array, p, MinPower(), ComplexWeights(), SOCP(), Clarabel.Optimizer)
```

## Related

**Alternative formulations:** [`LP`](@ref "LP"), [`QP`](@ref "QP")

**Objective requiring SOCP:** [`MaxDirectivity`](@ref "MaxDirectivity")

**Robust synthesis:** [`Robustness`](@ref)

**Background:** [Modulus Constraints per Formulation](@ref)
