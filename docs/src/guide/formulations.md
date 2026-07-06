# Formulations

The `AbstractFormulation` controls how the modulus constraint
$|AF(\hat{r}_p)| \leq t_p$ is enforced and which solver types are compatible. Pass it
as the fifth argument to [`synthesize`](@ref).

!!! note "Provisional guide"
    This page is a selection guide. The individual formulation pages contain syntax,
    examples, solver suggestions, and detailed notes.

## Choosing a Formulation

| Formulation | Modulus constraint | Objective type | Typical solver |
|---|---|---|---|
| [`LP`](@ref "LP") | Polyhedral, linear | Linear | HiGHS, Mosek |
| [`QP`](@ref "QP") | Polyhedral, linear | Linear or quadratic | HiGHS, Mosek |
| [`SOCP`](@ref "SOCP") | Exact second-order cone | Linear or quadratic | Clarabel, Mosek |

Use [`LP`](@ref "LP") as the default for linear objectives, especially when the array
factor is real. Use [`QP`](@ref "QP") when the objective has a quadratic term, such as
[`MinPower`](@ref "MinPower"), [`MinIntegratedPower`](@ref "MinIntegratedPower"), or
[`MinFieldError`](@ref "MinFieldError"). Use [`SOCP`](@ref "SOCP") when the exact
complex modulus constraint is important, or when an objective requires conic
constraints, such as [`MaxDirectivity`](@ref "MaxDirectivity").

For a real AF (`ConjugateSymmetricWeights` / `ProgressivePhaseAmplitude` +
`SymmetricArray`), the modulus constraint is $-t \leq AF \leq t$ in all three
formulations: exact and linear.

## Related

**Background:** [Modulus Constraints per Formulation](@ref)

**Excitations affecting the AF type:** [`ConjugateSymmetricWeights`](@ref "ConjugateSymmetricWeights"),
[`ProgressivePhaseAmplitude`](@ref "ProgressivePhaseAmplitude")
