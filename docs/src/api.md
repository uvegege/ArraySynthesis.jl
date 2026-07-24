# API Reference

## Synthesis

```@docs
synthesize
array_factor
SynthesisResult
IterativeSynthesisResult
MultiPatternResult
```

## Arrays

```@docs
ArrayGeometry
SymmetricArray
uniform_linear_array
linear_array
symmetric_linear_array
planar_array
triangular_array
hexagonal_array
circular_array
symmetric_planar_array
symmetric_triangular_array
symmetric_hexagonal_array
symmetric_circular_array
symmetrize
materialize
is_origin
is_representative
```

## Pattern

```@docs
AbstractDirection
ThetaDirection
UVDirection
Direction
Region
Pattern
Beam
ShapedBeam
NullPoint
Nulls
SideLobeRegion
ClosedInterval
..
direction
θ
uv
region
pattern
beam
shaped_beam
null
nulls
sidelobes
outside
join_regions
```

## Pattern Masks

```@docs
theta_ramp
u_ramp
v_ramp
csc_values
```

## Region Shapes

```@docs
RegionShape
Circle
Ellipse
Rectangle
Polygon
Moonlike
rhombus
triangle
visible_region
```

## Excitations

```@docs
AbstractExcitation
ComplexWeights
RealWeights
ConjugateSymmetricWeights
ProgressivePhaseAmplitude
```

## Formulations

```@docs
AbstractFormulation
LP
QP
SOCP
MILP
```

## Objectives

```@docs
AbstractObjective
DirectObjective
Feasible
MinSLL
MaxAF
MinL1
MinWeightedL1
MinPower
MinIntegratedPower
MinFieldError
```

## Iterative Methods

```@docs
SynthesisMethod
IterativeFloorSynthesis
IterativeReweightedL1
MultiPatternReweightedL1
```

## Robustness

```@docs
Tolerances
robust
```
