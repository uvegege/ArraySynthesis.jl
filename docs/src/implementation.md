# Implementation Notes

Short notes on the code.

## Common array-factor expressions

`array_factor_reim` builds common real/imaginary array-factor expressions for the
supported array and excitation types. Constraints and objectives then reuse those
expressions instead of maintaining separate matrix-form code paths for every
combination of array, excitation, formulation, and objective.

## Manual `AffExpr` construction

`polygon_bound!` builds each linear combination with `add_to_expression!`. This avoids the intermediate allocations from generic JuMP expression arithmetic in the innermost sidelobe-point / polygon-face loop.

## Solver choice for large problems

`HiGHS` and `Clarabel` are fine for small examples. Large planar problems or dense sampling grids are typically much faster with `Mosek`.
