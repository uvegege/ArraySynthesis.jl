# MaxDirectivity

Maximize directivity at a given direction relative to integrated power over a region.

## Syntax

```julia
MaxDirectivity(direction, region)
```

## Arguments

| Argument | Description |
|---|---|
| `direction` | Look direction (`ThetaDirection`, `UVDirection`, angle in radians, or `(u,v)` tuple). |
| `region` | Integration region for the denominator (total radiated power). |

## Description

Fixes $\mathrm{Re}(AF(\hat{r}_0)) = 1$ and minimizes the integrated power
$\sum_p |AF(\hat{r}_p)|^2$ over `region`, which is equivalent to maximizing
directivity $|AF(\hat{r}_0)|^2 / \int |AF|^2 \, d\Omega$.

Requires `SOCP`.

## Example

```julia
using ArraySynthesis; using ArraySynthesis: °, dB; using Clarabel

array  = uniform_linear_array(32, d = 0.5)
full_r = region(-90°..90°, 1°)
p      = pattern()   # no additional constraints

result = synthesize(array, p, MaxDirectivity(0°, full_r), ComplexWeights(), SOCP(), Clarabel.Optimizer)
```

## Related

**Related objectives:** [`MaxGain`](@ref "MaxGain"),
[`MinIntegratedPower`](@ref "MinIntegratedPower")

**Required formulation:** [`SOCP`](@ref "SOCP")
