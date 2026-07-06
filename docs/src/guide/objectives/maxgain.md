# MaxGain

Maximize the array gain at one or more directions.

## Syntax

```julia
MaxGain()
MaxGain(direction)
MaxGain(directions)
```

## Arguments

| Argument | Description |
|---|---|
| `direction` | A single direction (`ThetaDirection`, `UVDirection`, angle in radians, or `(u,v)` tuple). |
| `directions` | A vector of directions. |

If no argument is given, the beam directions defined in the `Pattern` are used.

## Description

Maximizes $\sum \mathrm{Re}(AF(\hat{r}_k))$ over the specified directions, subject to
the pattern constraints. Compatible with all formulations.

## Example

```julia
using ArraySynthesis; using ArraySynthesis: °, dB; using HiGHS

array = uniform_linear_array(20, d = 0.5)
p     = pattern(beam(-13°),
                sidelobes(region(-90°..(-4°), 1°), -40dB),
                sidelobes(region(15°..90°,    1°), -40dB))

# Direction taken from the beam in the pattern
result = synthesize(array, p, MaxGain(), ComplexWeights(), LP(), HiGHS.Optimizer)

# Explicit direction
result = synthesize(array, p, MaxGain(-13°), ComplexWeights(), LP(), HiGHS.Optimizer)
```

## Related

**Alternative objectives:** [`MinSLL`](@ref "MinSLL"), [`MaxDirectivity`](@ref "MaxDirectivity"),
[`Feasible`](@ref "Feasible")
