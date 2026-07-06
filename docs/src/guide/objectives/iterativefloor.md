# IterativeFloorSynthesis

Iterative QP for narrow-beam low-sidelobe synthesis.

## Syntax

```julia
IterativeFloorSynthesis(direction; sll = -30dB, max_iter = 15, tol = 1e-4, stable_iterations = 3)
```

## Arguments

| Argument | Description |
|---|---|
| `direction` | Look direction for the main beam. |
| `sll` | Target SLL in linear scale (default `-30dB`). |
| `max_iter` | Maximum number of iterations (default `15`). |
| `tol` | Convergence tolerance (default `1e-4`). |
| `stable_iterations` | Stop early if the active sidelobe set does not change for this many consecutive iterations (default `3`). |

## Description

Implements the Orchard–Elliott / Tseng–Griffiths iterative peak correction algorithm.
At each step, the highest sidelobe peak in the current solution is identified and a
linear floor constraint is added to the model. The process repeats until the SLL target
is reached or no further improvement occurs.

The sidelobe floor constraints grow with each iteration, increasing model size.
Requires `QP`.

## Example

```julia
using ArraySynthesis; using ArraySynthesis: °, dB; using HiGHS

array = uniform_linear_array(20, d = 0.5)
p     = pattern(sidelobes(region(-90°..34°, 0.5°), -40dB),
                sidelobes(region(74°..90°,  0.5°), -40dB))
obj   = IterativeFloorSynthesis(ThetaDirection(50°); sll = -40dB, max_iter = 8)

result = synthesize(array, p, obj, ComplexWeights(), QP(), HiGHS.Optimizer)
```

## Related

**Related objectives:** [`IterativePatternLeastSquares`](@ref "IterativePatternLeastSquares"),
[`MinSLL`](@ref "MinSLL")

**Required formulation:** [`QP`](@ref "QP")
