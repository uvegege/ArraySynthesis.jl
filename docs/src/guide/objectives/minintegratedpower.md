# MinIntegratedPower

Minimize integrated power over one or more regions.

## Syntax

```julia
MinIntegratedPower(region)
MinIntegratedPower(regions)
```

## Arguments

| Argument | Description |
|---|---|
| `region` / `regions` | A `Region` or a vector of `Region` objects over which power is integrated. |

## Description

Minimizes $\sum_p |AF(\hat{r}_p)|^2$ over the discretized directions in the specified
regions, subject to the pattern constraints. The objective is quadratic; requires
`QP` or `SOCP`.

Useful for interference suppression over a sector or for minimizing the energy radiated
in an unwanted direction, while maintaining gain and sidelobe constraints elsewhere.

## Example

```julia
using ArraySynthesis; using ArraySynthesis: °, dB; using HiGHS

array      = uniform_linear_array(32, d = 0.5)
interf_r   = region(-30°..(-10°), 1°)
p          = pattern(beam(0°), sidelobes(region(-90°..(-35°), 1°), -20dB),
                                sidelobes(region(5°..90°,     1°), -20dB))

result = synthesize(array, p, MinIntegratedPower(interf_r), ComplexWeights(), QP(), HiGHS.Optimizer)
```

## Related

**Alternative objectives:** [`MinPower`](@ref "MinPower"), [`MinSLL`](@ref "MinSLL")

**Compatible formulations:** [`QP`](@ref "QP"), [`SOCP`](@ref "SOCP")
