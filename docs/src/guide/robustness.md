# Robustness

Robustness lets the model reserve margin against small pointing, channel, and
position errors. It is enabled with the `robustness` keyword in [`synthesize`](@ref):

```julia
rob = robust(
    pointing_accuracy = 0.1°,
    phase_tolerance = 0.5°,
    amplitude_tolerance = 0.05dB,
    position_tolerance = (xy = 0.0005, z = 0.0),
)

result = synthesize(array, p, MinPower(), ComplexWeights(), SOCP(), Clarabel.Optimizer;
                    robustness = rob)
```

The nominal bound

```math
|AF(p)| \leq U(p)
```

is replaced by

```math
|AF(p)| + \delta(p)\|w\|_2 \leq U(p).
```

Internally this is built as one shared cone for ``\|w\|_2`` and an affine correction
``\delta(p)t`` in each affected bound. This keeps the implementation small and avoids
adding one cone per sampled point.

## What Is Robustified

Robustness currently applies to:

- `sidelobes`
- `null`
- `shaped_beam`
- [`MinSLL`](@ref "MinSLL")

Beam equalities are kept nominal and act as normalization constraints.
The current implementation is wired for direct single-shot solves. Iterative
strategies could also pass the same robust bounds to their inner problems, but
that path is not connected yet.

## Tolerances

`robust` accepts physical tolerances:

- `pointing_accuracy`: angular error, in radians for `ThetaDirection`; for
  `UVDirection` it is interpreted as an offset in `(u, v)` coordinates.
- `phase_tolerance`: per-channel phase error, in radians.
- `amplitude_tolerance`: per-channel amplitude error. Values such as `0.05dB`
  are converted to a linear relative tolerance.
- `position_tolerance`: position error in wavelengths. A scalar means isotropic
  3D error; `(xy = ..., z = ...)` separates in-plane and out-of-plane errors.

The default combination is `combine = :quadrature`; use `combine = :sum` for a
more conservative direct sum.

## Formulation

This implementation is SOCP-only because the correction uses ``\|w\|_2``. Calling
`synthesize(...; robustness = rob)` with `LP()` or `QP()` raises an error.

Deep nulls and tight sidelobe masks can become infeasible when the uncertainty
margin is larger than the available bound. That is expected: robust synthesis is
trading nominal performance for certified margin.

## Example

See `examples/SOCP/robust_sidelobes.jl`.
