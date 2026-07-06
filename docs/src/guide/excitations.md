# Excitations

The `AbstractExcitation` defines the parametrization of the array weights and,
consequently, the structure of the array factor passed to the optimizer. Pass it as
the fourth argument to [`synthesize`](@ref).

!!! note "Provisional guide"
    This page is a selection guide. The individual excitation pages contain the
    syntax, examples, and detailed notes.

## Choosing an Excitation

| Excitation | Array type | Variables | Best suited for |
|---|---|---|---|
| [`ComplexWeights`](@ref "ComplexWeights") | `ArrayGeometry` | Complex weights | General unconstrained synthesis |
| [`RealWeights`](@ref "RealWeights") | `ArrayGeometry` | Real weights | Fixed-phase or amplitude-only feeding |
| [`ConjugateSymmetricWeights`](@ref "ConjugateSymmetricWeights") | `SymmetricArray` | Complex half-array weights | Real array factor with maximum flexibility |
| [`ProgressivePhaseAmplitude`](@ref "ProgressivePhaseAmplitude") | Any array | Real amplitudes with fixed phase gradient | Steered beams with known phase progression |

## Real Array Factors

[`ConjugateSymmetricWeights`](@ref "ConjugateSymmetricWeights") is the main option for
centrosymmetric arrays. It enforces conjugate-symmetric weights and makes the array
factor real by construction. On a `SymmetricArray`,
[`ProgressivePhaseAmplitude`](@ref "ProgressivePhaseAmplitude") also gives a real array
factor while fixing the phase progression toward a steering direction.

A real array factor is important because the modulus constraint becomes a pair of
linear inequalities. See [Symmetric Arrays and Real Array Factor](@ref) for the
derivation.

## Related

**Background:** [Excitation Types](@ref), [Symmetric Arrays and Real Array Factor](@ref)

**Formulation impact:** [`LP`](@ref "LP"), [`QP`](@ref "QP"), [`SOCP`](@ref "SOCP")
