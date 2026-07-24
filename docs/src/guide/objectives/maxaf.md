# MaxAF

Maximize the real array-factor response at one or more directions.

## Syntax

```julia
MaxAF()
MaxAF(direction)
MaxAF(directions)
```

## Arguments

| Argument | Description |
|---|---|
| `direction` | A single direction (`ThetaDirection`, `UVDirection`, angle in radians, or `(u,v)` tuple). |
| `directions` | A vector of directions. |

If no argument is given, the beam directions defined in the `Pattern` are used.

## Description

Maximizes $\sum \mathrm{Re}(AF(\hat{r}_k))$ over the specified directions, subject to the pattern constraints. Compatible with all formulations.

## Example

```julia
using ArraySynthesis; using ArraySynthesis: °, dB; using HiGHS

array = uniform_linear_array(20, d = 0.5)
p = pattern(beam(-13°),
            sidelobes(region(-90°..(-4°), 1°), -40dB),
            sidelobes(region(15°..90°,    1°), -40dB))

# Direction taken from the beam in the pattern
result = synthesize(array, p, MaxAF(), ComplexWeights(), LP(), HiGHS.Optimizer)

# Explicit direction
result = synthesize(array, p, MaxAF(-13°), ComplexWeights(), LP(), HiGHS.Optimizer)
```
