# MinPower

Minimize total excitation power subject to pattern constraints.

## Syntax

```julia
MinPower()
```

## Description

Minimizes $\sum_n |w_n|^2$, the total radiated power, subject to the beam, sidelobe,
and null constraints defined in the `Pattern`. The objective is quadratic; requires
`QP` or `SOCP`.

A common use case is minimum-power beamforming: fix the gain at a look direction
(`beam`) and minimize power subject to sidelobe and null constraints.

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

**Power objectives:** [`MinIntegratedPower`](@ref "MinIntegratedPower"),
[`MaxDirectivity`](@ref "MaxDirectivity")

**Compatible formulations:** [`QP`](@ref "QP"), [`SOCP`](@ref "SOCP")
