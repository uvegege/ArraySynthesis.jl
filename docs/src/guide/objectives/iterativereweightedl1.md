# IterativeReweightedL1

Sparse array synthesis via iterative reweighted L1 minimization.

## Syntax

```julia
IterativeReweightedL1(; max_iter = 15, d0 = 1e-3, stable_iterations = 3)
```

## Arguments

| Argument | Description |
|---|---|
| `max_iter` | Maximum number of iterations (default `15`). |
| `d0` | Initial regularization to avoid division by zero when computing weights (default `1e-3`). |
| `stable_iterations` | Stop early if the active element set does not change for this many consecutive iterations (default `3`). |

## Description

Promotes sparse arrays using iterative reweighted L1 minimization. At each iteration
$k$, elements with small excitation amplitude $|w_n^{(k)}|$ receive a higher penalty
weight $\alpha_n^{(k+1)} \propto 1 / (|w_n^{(k)}| + d_0)$ in the next L1 minimization.
This biases the solution toward fewer active elements without requiring integer
variables.

The number of active elements is not set explicitly; it emerges from the pattern
constraints and the array geometry. Compatible with `LP` and `SOCP`.

## Literature

The reweighting idea is the Candès–Wakin–Boyd sparse recovery heuristic. In antenna
array synthesis, this package's single-pattern sparse strategy is closest to the
sequential convex optimization approach in:

- B. Fuchs, "Synthesis of Sparse Arrays With Focused or Shaped Beampattern via
  Sequential Convex Optimizations," *IEEE Transactions on Antennas and Propagation*,
  vol. 60, no. 7, pp. 3499-3503, Jul. 2012.

## Example

```julia
using ArraySynthesis; using ArraySynthesis: °, dB; using Clarabel

array  = triangular_array(17, 17, d = 0.5)
beam_r = region(Circle(0.27, (0.0, 0.0)), step = 2°)
sl_r   = visible_region(Circle(0.27, (0.0, 0.0)); step = 4°, bandpass = 0.2)
p      = pattern(shaped_beam(beam_r, 1.0, ripple = -1.5dB), sidelobes(sl_r, -25dB))
obj    = IterativeReweightedL1(max_iter = 15)

result = synthesize(array, p, obj, ComplexWeights(), SOCP(), Clarabel.Optimizer)
```

Check how many elements are active after synthesis:

```julia
active = abs.(result.weights) .> 1e-3
println("Active: $(sum(active)) / $(length(active))")
```

## Related

**Sparse objectives:** [`MultiPatternReweightedL1`](@ref "MultiPatternReweightedL1"),
[`IterativePatternLeastSquares`](@ref "IterativePatternLeastSquares")
