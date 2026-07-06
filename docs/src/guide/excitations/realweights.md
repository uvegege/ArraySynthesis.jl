# RealWeights

Real-valued excitations, one per element.

## Syntax

```julia
RealWeights()
```

## Description

Each element has a real weight $w_n \in \mathbb{R}$. The array factor is still complex
in general (because the steering matrix $A$ is complex), but the number of optimization
variables is halved compared to `ComplexWeights`.

Appropriate when the feeding network physically constrains all excitations to be
real-valued (e.g. uniform-phase corporate feed, Taylor distribution starting point).

Not to be confused with `ConjugateSymmetricWeights`, which enforces a symmetry
condition that makes the **array factor** real. `RealWeights` only constrains
the weights, not the AF.

- Use with `ArrayGeometry`.
- Compatible with all formulations.

## Related

**Alternative excitations:** [`ComplexWeights`](@ref "ComplexWeights"),
[`ConjugateSymmetricWeights`](@ref "ConjugateSymmetricWeights")
