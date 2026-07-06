# ProgressivePhaseAmplitude

Real amplitudes with a fixed progressive phase gradient toward a steering direction.

## Syntax

```julia
ProgressivePhaseAmplitude()
ProgressivePhaseAmplitude(β)
```

## Arguments

| Argument | Description |
|---|---|
| `β` | Steering direction. Accepts a `ThetaDirection`, a `UVDirection`, a single angle (radians), a `(u, v)` tuple, or `nothing`. If omitted or `nothing`, the direction is extracted from the main beam in the `Pattern`. |

## Description

Parametrizes the excitation as $w_n = a_n e^{j\boldsymbol{\beta}\cdot\mathbf{p}_n}$,
where $\boldsymbol{\beta}$ is a fixed phase gradient and $a_n \in \mathbb{R}$ are the
optimization variables. The steering matrix is evaluated relative to $\boldsymbol{\beta}$,
so constraints at the look direction are exact and linear regardless of formulation.

On a `SymmetricArray`, the array factor is real (same as `ConjugateSymmetricWeights`).
On a general `ArrayGeometry`, the array factor is complex.

## Comparison with `ConjugateSymmetricWeights`

| | `ConjugateSymmetricWeights` | `ProgressivePhaseAmplitude` |
|---|---|---|
| Variables | $\mathrm{Re}(w_n)$, $\mathrm{Im}(w_n)$ (complex, half array) | $a_n$ (real amplitudes) |
| Phase | Free per element pair | Fixed by $\boldsymbol{\beta}$ |
| Real AF | Yes (on `SymmetricArray`) | Yes (on `SymmetricArray`) |
| Best for | Maximum pattern flexibility | Known look direction, amplitude-only control |

## Example

```julia
using ArraySynthesis; using ArraySynthesis: °, dB; using HiGHS

array  = symmetric_linear_array(32, d = 0.5)
beam_r = region(12.5°..37.5°, 1°)
sll_r  = join_regions(region(-90°..6.5°, 1°), region(43.5°..90°, 1°))
p      = pattern(shaped_beam(beam_r, 1.0, ripple = -0.6dB))

# Explicit steering direction
result = synthesize(array, p, MinSLL(sll_r), ProgressivePhaseAmplitude(deg2rad(25)), LP(), HiGHS.Optimizer)

# Extract steering direction from pattern automatically
result = synthesize(array, p, MinSLL(sll_r), ProgressivePhaseAmplitude(), LP(), HiGHS.Optimizer)
```

## Related

**Alternative excitations:** [`ConjugateSymmetricWeights`](@ref "ConjugateSymmetricWeights"),
[`ComplexWeights`](@ref "ComplexWeights")

**Background:** [Excitation Types](@ref)
