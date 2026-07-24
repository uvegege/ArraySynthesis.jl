# Formulations

The `AbstractFormulation` controls how the modulus constraint
$|AF(\hat{r}_p)| \leq t_p$ is enforced and which solver types are compatible. Pass it
as the fifth argument to [`synthesize`](@ref).

Use this page to choose the optimization formulation. For exact constructor
signatures, see the [API Reference](@ref).

## Choosing a Formulation

| Formulation | Modulus constraint | Objective type | Typical solver |
|---|---|---|---|
| [`LP`](@ref "LP") | Linear | Linear | HiGHS |
| [`QP`](@ref "QP") | Linear | Linear or quadratic | HiGHS |
| [`SOCP`](@ref "SOCP") | Exact second-order cone | Linear or quadratic | Clarabel |

Use [`LP`](@ref "LP") as the default for linear objectives, especially when the array
factor is real. Use [`QP`](@ref "QP") when the objective has a quadratic term, such as
[`MinPower`](@ref "MinPower"), [`MinIntegratedPower`](@ref "MinIntegratedPower"), or
[`MinFieldError`](@ref "MinFieldError"). Use [`SOCP`](@ref "SOCP") when exact
complex magnitude constraints or robust synthesis margins are important.

For a real AF (`ConjugateSymmetricWeights` / `ProgressivePhaseAmplitude` +
`SymmetricArray`), the modulus constraint is $-t \leq AF \leq t$ in all three
formulations: exact and linear.

## Formulation Notes

[`LP`](@ref "LP") is the default for linear objectives. For complex array factors, magnitude constraints are approximated with a polygon controlled by `polygon_faces`.

[`QP`](@ref "QP") uses the same polygonal magnitude constraints as `LP`, but allows quadratic objectives such as [`MinPower`](@ref "MinPower"), [`MinIntegratedPower`](@ref "MinIntegratedPower"), and [`MinFieldError`](@ref "MinFieldError").

[`SOCP`](@ref "SOCP") represents complex magnitude bounds exactly with second-order cones. Use it for exact complex masks and robust synthesis.
