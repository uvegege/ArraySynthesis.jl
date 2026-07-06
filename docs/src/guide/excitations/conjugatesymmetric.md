# ConjugateSymmetricWeights

Conjugate-symmetric excitations for centrosymmetric arrays, producing a real array factor.

## Syntax

```julia
ConjugateSymmetricWeights()
```

## Description

For a centrosymmetric array (positions in symmetric pairs $\pm\mathbf{p}_n$), enforcing
$w_{-n} = w_n^*$ makes the array factor **always real**:

```math
AF(\hat{r}) = \sum_n 2\cos(2\pi\hat{r}\cdot\mathbf{p}_n)\,\mathrm{Re}(w_n)
            - 2\sin(2\pi\hat{r}\cdot\mathbf{p}_n)\,\mathrm{Im}(w_n)
```

Only the $\lceil N/2 \rceil$ representative weights (stored in `SymmetricArray`) are
optimization variables.

**Advantages over `ComplexWeights`:**
- The modulus constraint $|AF| \leq t$ reduces to $-t \leq AF \leq t$: exact and linear
  in all formulations, with no polygon approximation error.
- Fewer optimization variables.

Requires `SymmetricArray` (`symmetric_linear_array`, `symmetric_planar_array`, etc.).

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

**Alternative excitations:** [`ProgressivePhaseAmplitude`](@ref "ProgressivePhaseAmplitude"),
[`RealWeights`](@ref "RealWeights")

**Background:** [Symmetric Arrays and Real Array Factor](@ref)
