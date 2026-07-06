# MultiPatternReweightedL1

Sparse array synthesis for multiple patterns sharing the same active element support.

## Syntax

```julia
MultiPatternReweightedL1(; max_iter = 15, d0 = 1e-3, stable_iterations = 3)
```

## Arguments

| Argument | Description |
|---|---|
| `max_iter` | Maximum number of iterations (default `15`). |
| `d0` | Initial regularization (default `1e-3`). |
| `stable_iterations` | Stop early if the active set is stable (default `3`). |

## Description

Extension of `IterativeReweightedL1` to reconfigurable arrays that must satisfy $K$
independent patterns with a single physical aperture (same set of active elements).
A shared L1 auxiliary variable couples all patterns in the sparsity objective; each
pattern has independent weight variables optimized separately.

Pass a **vector** of `Pattern` objects as the second argument to `synthesize`.
The result is a `MultiPatternResult` with one weight vector per pattern.

## Example

Two-pattern reconfigurable planar array — broadside circular beam and offset rhombus beam:

```julia
using ArraySynthesis; using ArraySynthesis: °, dB; using HiGHS

array = planar_array(14, 14, dx = 0.5, dy = 0.5)

p1 = pattern(
    shaped_beam(region(Circle(0.2, (0.0, 0.0)), step = 2°), 1.0, ripple = -1dB),
    sidelobes(visible_region(Circle(0.4, (0.0, 0.0)); step = 4°, bandpass = 0.0), -25dB))

p2 = pattern(
    shaped_beam(region(rhombus((0.2, 0.2), 0.2), step = 2°), 1.0, ripple = -1dB),
    sidelobes(visible_region(rhombus((0.2, 0.2), 0.4); step = 4°, bandpass = 0.0), -24dB))

obj    = MultiPatternReweightedL1(max_iter = 15)
result = synthesize(array, [p1, p2], obj, ComplexWeights(), LP(), HiGHS.Optimizer)

# Active elements: used in any of the K patterns
max_activity = [maximum(abs(result.weights[k][n]) for k in eachindex(result.weights))
                for n in eachindex(result.weights[1])]
active = max_activity .> 1e-5
```

## Related

**Sparse objectives:** [`IterativeReweightedL1`](@ref "IterativeReweightedL1")
