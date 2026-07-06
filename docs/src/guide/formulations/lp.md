# LP

Linear programming formulation with polyhedral modulus constraint.

## Syntax

```julia
LP()
LP(polygon_faces = F)
```

## Arguments

| Argument | Description |
|---|---|
| `polygon_faces` | Number of polygon faces for the modulus approximation (default `8`). |

## Description

Approximates $|AF(\hat{r}_p)| \leq t_p$ by a regular $F$-faced inner polygon:

```math
\mathrm{Re}(e^{j\theta_\ell} AF(\hat{r}_p)) \leq t_p \cos(\pi/F),
\quad \theta_\ell = \frac{2\pi\ell}{F}, \quad \ell = 0,\ldots,F-1
```

The resulting model is a linear program solvable by any LP solver.

The approximation is conservative (inner polygon); the relative tightening is
$1 - \cos(\pi/F)$. With the default $F = 8$ this is ~7.6%.
When the AF is real (`ConjugateSymmetricWeights` or `ProgressivePhaseAmplitude` on a
`SymmetricArray`), the constraint reduces to $-t_p \leq AF \leq t_p$, which is exact
and linear regardless of `polygon_faces`.

Suggested solvers: `HiGHS.Optimizer`, `Mosek.Optimizer`

## Example

```julia
using ArraySynthesis; using ArraySynthesis: °, dB; using HiGHS

beam_r = region(12.5°..37.5°, 1°)
sll_r  = join_regions(region(-90°..6.5°, 1°), region(43.5°..90°, 1°))
p      = pattern(shaped_beam(beam_r, 1.0, ripple = -0.6dB))
array  = symmetric_linear_array(32, d = 0.5)

result = synthesize(array, p, MinSLL(sll_r), ConjugateSymmetricWeights(), LP(), HiGHS.Optimizer)
```

## Related

**Alternative formulations:** [`QP`](@ref "QP"), [`SOCP`](@ref "SOCP")

**Background:** [Modulus Constraints per Formulation](@ref)
