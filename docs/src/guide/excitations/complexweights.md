# ComplexWeights

Fully free complex excitations, one per element.

## Syntax

```julia
ComplexWeights()
```

## Description

Each element has an independent complex weight $w_n \in \mathbb{C}$. No structure is
imposed on the weights; this is the most general excitation model.

The array factor is complex and the modulus constraint $|AF(\hat{r}_p)| \leq t_p$ is
approximated by a polygon (LP, QP) or imposed exactly as a second-order cone (SOCP).

- Use with `ArrayGeometry` or `SymmetricArray` (via `materialize`).
- Compatible with all formulations and all objectives.

## Example

```julia
using ArraySynthesis; using ArraySynthesis: °, dB; using HiGHS

array = uniform_linear_array(20, d = 0.5)
p     = pattern(beam(-13°),
                sidelobes(region(-90°..(-4°), 1°), -40dB),
                sidelobes(region(15°..90°,    1°), -40dB))

result = synthesize(array, p, MaxGain(), ComplexWeights(), LP(), HiGHS.Optimizer)
```

## Related

**Alternative excitations:** [`RealWeights`](@ref "RealWeights"),
[`ConjugateSymmetricWeights`](@ref "ConjugateSymmetricWeights"),
[`ProgressivePhaseAmplitude`](@ref "ProgressivePhaseAmplitude")
