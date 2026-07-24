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

Compatible with linear and quadratic objectives. Typically slower than `LP`/`QP` on large models but more accurate for complex AF.

Also required when using [`robust`](@ref), because robust bounds use one shared norm cone for the excitation vector.

Suggested solvers: `Clarabel.Optimizer`, `Mosek.Optimizer`

## Example

```julia
using ArraySynthesis; using ArraySynthesis: °, dB; using Clarabel

array = uniform_linear_array(32, d = 0.5)
p = pattern(beam(0°),
            sidelobes(region(-90°..(-5°), 1°), -20dB),
            sidelobes(region(5°..90°,     1°), -20dB),
            null(-40°, level = -55dB))

result = synthesize(array, p, MinPower(), ComplexWeights(), SOCP(), Clarabel.Optimizer)
```