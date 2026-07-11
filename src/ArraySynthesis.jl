module ArraySynthesis

    using LinearAlgebra
    using JuMP

    include("./pattern.jl")
    include("./pattern_mask.jl")
    include("./pattern_shapes.jl")
    include("./array.jl")
    include("./coefficients.jl")
    include("./steering.jl")
    include("./variables.jl")
    include("./robustness.jl")
    include("./constraints.jl")
    include("./objective.jl")
    include("./strategy.jl")
    include("./synthesize.jl")


    export synthesize, array_factor

    export ArrayGeometry, SymmetricArray
    export uniform_linear_array, linear_array, symmetric_linear_array
    export planar_array, triangular_array, hexagonal_array, circular_array
    export symmetric_planar_array, symmetric_triangular_array
    export symmetric_hexagonal_array, symmetric_circular_array
    export symmetrize, materialize, is_origin, is_representative

    export AbstractDirection, ThetaDirection, UVDirection, Direction
    export Region, Pattern, Beam, ShapedBeam, NullPoint, Nulls, SideLobeRegion
    export ClosedInterval, ..
    export uv, θ, direction
    export region, pattern, beam, shaped_beam, sidelobes, null, nulls
    export outside, join_regions

    export theta_ramp, u_ramp, v_ramp, csc_values

    export RegionShape, Circle, Ellipse, Rectangle, Polygon, Moonlike
    export rhombus, triangle, visible_region

    export AbstractExcitation
    export ComplexWeights, RealWeights, ConjugateSymmetricWeights, ProgressivePhaseAmplitude

    export AbstractFormulation, LP, QP, SOCP, MILP

    export AbstractObjective, DirectObjective
    export Feasible, MaxGain, MinPower, MinSLL, MinIntegratedPower
    export MinL1, MinWeightedL1, MinFieldError, MaxDirectivity

    export SynthesisMethod
    export IterativeFloorSynthesis, IterativePatternLeastSquares, IterativeReweightedL1
    export MultiPatternReweightedL1

    export SynthesisResult, IterativeSynthesisResult, MultiPatternResult

    export Tolerances, robust

end
