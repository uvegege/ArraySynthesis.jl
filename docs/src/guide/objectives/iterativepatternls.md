# IterativePatternLeastSquares

Iterative least-squares for shaped beam synthesis.

## Syntax

```julia
IterativePatternLeastSquares(region, target; max_iter = 15)
```

## Arguments

| Argument | Description |
|---|---|
| `region` | `Region` defining the shaped beam area. |
| `target` | Target amplitude values: a vector of the same length as `region.points`, or a scalar. |
| `max_iter` | Maximum number of iterations (default `15`). |

## Description

Based on Fuchs (2010). At each iteration, the phase reference is updated from the
previous solution, so the optimizer fits a complex target without requiring the phase
to be fixed in advance. The method converges to a shaped beam that matches the target
amplitude while satisfying the sidelobe and null constraints from the `Pattern`.

Unlike `MinFieldError`, the target phase is not fixed: the algorithm finds the phase
that best matches the target amplitude, making it suitable for arbitrary arrays where
a natural phase reference does not exist.

Compatible with `LP` and `SOCP`.

## Example

Cosine-shaped beam over $[-25°, 25°]$ with -25 dB sidelobes:

```julia
using ArraySynthesis; using ArraySynthesis: °, dB; using Clarabel

array  = uniform_linear_array(28, d = 0.5)
beam_r = region(-25°..25°, 1°)
sll_r  = join_regions(region(-90°..(-32°), 1°), region(32°..90°, 1°))
target = [0.55 + 0.45 * cos(π * p.θ / 25°) for p in beam_r.points]
p      = pattern(sidelobes(sll_r, -25dB))
obj    = IterativePatternLeastSquares(beam_r, target; max_iter = 12)

result = synthesize(array, p, obj, ComplexWeights(), SOCP(), Clarabel.Optimizer)
```

## Related

**Field-fitting objectives:** [`MinFieldError`](@ref "MinFieldError"),
[`IterativeFloorSynthesis`](@ref "IterativeFloorSynthesis"),
[`MinSLL`](@ref "MinSLL")
