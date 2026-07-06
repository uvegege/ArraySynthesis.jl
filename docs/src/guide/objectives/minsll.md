# MinSLL

Minimize the peak sidelobe level over one or more regions.

## Syntax

```julia
MinSLL(region; lower_bound = 0.0, upper_bound = -20dB)
MinSLL(regions; lower_bound = 0.0, upper_bound = -20dB)
```

## Arguments

| Argument | Description |
|---|---|
| `region` / `regions` | A `Region` or a vector of `Region` objects defining the sidelobe area. |
| `lower_bound` | Lower bound on the SLL auxiliary variable (default `0.0`). |
| `upper_bound` | Upper bound on the SLL auxiliary variable (default `-20 dB`). Tighten if infeasibility occurs at the default. |

## Description

Introduces one auxiliary variable $t_i$ per region and minimizes $\sum_i t_i$ subject
to $|AF(\hat{r}_p)| \leq t_i$ for each discretized direction $\hat{r}_p$ in region $i$.
Compatible with all formulations.

When the AF is real (`ConjugateSymmetricWeights` or `ProgressivePhaseAmplitude` on a
`SymmetricArray`), the sidelobe constraint reduces to $-t_i \leq AF(\hat{r}_p) \leq t_i$,
which is exact and linear with no polygon approximation error.

## Example

Flat-top shaped beam with minimized sidelobe level:

```julia
using ArraySynthesis; using ArraySynthesis: °, dB; using HiGHS

beam_r = region(12.5°..37.5°, 1°)
sll_r  = join_regions(region(-90°..6.5°, 1°), region(43.5°..90°, 1°))
p      = pattern(shaped_beam(beam_r, 1.0, ripple = -0.6dB))
array  = symmetric_linear_array(32, d = 0.5)

result = synthesize(array, p, MinSLL(sll_r), ConjugateSymmetricWeights(), LP(), HiGHS.Optimizer)
```

Multiple independent sidelobe regions:

```julia
obj = MinSLL([sll_r1, sll_r2, sll_r3])
```

## Related

**Alternative objectives:** [`MaxGain`](@ref "MaxGain"),
[`MinIntegratedPower`](@ref "MinIntegratedPower")

**Useful excitations:**
[`ConjugateSymmetricWeights`](@ref "ConjugateSymmetricWeights"),
[`ProgressivePhaseAmplitude`](@ref "ProgressivePhaseAmplitude")
