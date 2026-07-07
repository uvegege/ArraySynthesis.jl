# Implementation Notes

A few notes on how the model-building code is structured, mostly relevant if you're
reading the source or profiling a large synthesis problem.

## Dispatch over matrix branching

`array_factor_reim` is not one generic function that branches on excitation/array type
at runtime. It is a small set of methods dispatched on `AbstractExcitation` (and on
`ArrayGeometry` vs `SymmetricArray`), each assembling the real/imaginary array factor
differently (complex, real-only, amplitude-only, conjugate-symmetric). Julia resolves
the right method per call site, so there is no generic matrix path with `if`/`else`
branches on excitation type to maintain or pay for at every point evaluation.

## Manual `AffExpr` construction

The polygonal modulus bound (`polygon_bound!` in `constraints.jl`) builds the linear
combination of the real/imaginary `AffExpr`s term-by-term with `add_to_expression!`
instead of writing the natural `cos(θ)*re - sin(θ)*im` and letting JuMP's operator
overloading assemble it. This constraint is generated per sidelobe point per polygon
face, so for large regions or planar arrays it dominates build time; the manual form
avoids the intermediate allocations JuMP's generic expression arithmetic would produce.

## Solver choice for large problems

`HiGHS` and `Clarabel` are fine for the examples in this repo, but for large arrays
(planar, many sidelobe points, or MILP thinning) a commercial solver such as `Mosek`
is noticeably faster and more robust. See the `Formulations` guide for which solver
each formulation supports.
