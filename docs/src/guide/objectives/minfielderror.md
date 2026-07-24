# MinFieldError

Minimize the L2 error between the array factor and a complex reference pattern.

## Syntax

```julia
MinFieldError(region, reference)
```

## Arguments

| Argument | Description |
|---|---|
| `region` | A `Region` defining the evaluation points. |
| `reference` | Target complex values. A vector of the same length as `region.points`, a scalar, or a real number. |

## Description

Minimizes $\sum_p |AF(\hat{r}_p) - f_p|^2$ over the points in `region`, where $f_p$ are the complex reference values. The objective is quadratic; requires `QP` or `SOCP`.

The reference phase is fixed. Use `shaped_beam` when the specification is an amplitude mask with ripple rather than a complex reference.

## Example

```julia
using ArraySynthesis; using ArraySynthesis: °, dB; using HiGHS

array  = uniform_linear_array(28, d = 0.5)
beam_r = region(-20°..20°, 1°)
target = [cos(π * p.θ / 20°) for p in beam_r.points]   # real-valued cosine target
p = pattern(sidelobes(region(-90°..(-25°), 1°), -25dB),
            sidelobes(region(25°..90°,     1°), -25dB))

result = synthesize(array, p, MinFieldError(beam_r, target), ComplexWeights(), QP(), HiGHS.Optimizer)
```
