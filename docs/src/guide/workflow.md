# User Guide

ArraySynthesis problems are built from five pieces:

```julia
result = synthesize(array, pattern, objective, excitation, formulation, solver)
```

Each piece answers a different question.

| Piece | Question | Common choices |
|---|---|---|
| `array` | Where are the elements? | `uniform_linear_array`, `planar_array`, `symmetric_linear_array` |
| `pattern` | What should the array factor satisfy? | `pattern(beam(...), sidelobes(...), nulls(...))` |
| `objective` | What should be optimized? | `MinSLL`, `MaxAF`, `MinPower`, `MinL1` |
| `excitation` | How are weights parametrized? | `ComplexWeights`, `ConjugateSymmetricWeights` |
| `formulation` | How are magnitude constraints represented? | `LP`, `QP`, `SOCP` |

## Arrays

Array coordinates are expressed in wavelengths and stored as a `3 x N` matrix with one column per element. Convenience constructors cover the most common layouts:

```julia
array = uniform_linear_array(32, d = 0.5)
array = planar_array(16, 16, dx = 0.5, dy = 0.5)
array = circular_array([0.0, 1.0, 2.0], [1, 8, 16])
```

Use `SymmetricArray` constructors when the geometry is centrosymmetric and the
excitation model can exploit that symmetry:

```julia
array = symmetric_linear_array(32, d = 0.5)
coef = ConjugateSymmetricWeights()
```

This stores only one representative of each symmetric pair and gives a real array factor for conjugate-symmetric weights.

## Patterns

The pattern is the constraint specification. It is intentionally written as a small DSL so that the problem statement stays close to the antenna specification:

```julia
p = pattern(
    beam(0°),
    nulls([-35°, 35°]),
    sidelobes(region(-90°..-10°, 1°), -25dB),
    sidelobes(region(10°..90°, 1°), -25dB),
)
```

See [Pattern and Regions](@ref) for the detailed guide, including planar and conformal-array regions with [`visible_region`](@ref).

## Solving

Linear objectives and real array factors usually work well with `LP()` and HiGHS. Quadratic objectives need `QP()`. Exact complex magnitude constraints and robust synthesis use `SOCP()`.

```julia
using ArraySynthesis
using ArraySynthesis: °, dB
using HiGHS

array = symmetric_linear_array(32, d = 0.5)
sll = join_regions(region(-90°..6.5°, 1°), region(43.5°..90°, 1°))
p = pattern(
    shaped_beam(region(12.5°..37.5°, 1°), 1.0, ripple = -0.6dB),
)

result = synthesize(array, p, MinSLL(sll),
                    ConjugateSymmetricWeights(), LP(), HiGHS.Optimizer)
```

The returned result contains the synthesized complex weights:

```julia
w = result.weights
```

Use [`array_factor`](@ref) to evaluate the resulting pattern on any sampled grid.
