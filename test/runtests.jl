using ArraySynthesis
using ArraySynthesis: °, dB
using Test

@testset "array constructors" begin
    ula = uniform_linear_array(4, d = 0.5)
    @test ula isa ArrayGeometry
    @test size(ula.positions) == (3, 4)
    @test ula.dim == 1
    @test ula.positions[1, :] ≈ [-0.75, -0.25, 0.25, 0.75]

    lin = linear_array([-1.0, 0.0, 1.0], axis = :y)
    @test lin.positions[2, :] == [-1.0, 0.0, 1.0]

    planar = planar_array(2, 3, dx = 0.5, dy = 0.25)
    @test planar isa ArrayGeometry
    @test size(planar.positions) == (3, 6)
    @test planar.dim == 2

    triangular = triangular_array(2, 2)
    hexagonal = hexagonal_array(1)
    circular = circular_array([0.0, 1.0], [1, 4])
    @test size(triangular.positions, 2) == 4
    @test size(hexagonal.positions, 2) == 7
    @test size(circular.positions, 2) == 5
end

@testset "symmetric arrays" begin
    sym = symmetric_linear_array(5, d = 0.5)
    @test sym isa SymmetricArray
    @test size(sym.positions, 2) == 3
    @test any(is_origin, eachcol(sym.positions))

    full = materialize(sym)
    @test full isa ArrayGeometry
    @test size(full.positions, 2) == 5

    p, dim = symmetrize(full.positions)
    @test dim == 1
    @test size(p, 2) == size(sym.positions, 2)
end

@testset "directions and one-dimensional patterns" begin
    @test direction(0.0) isa ThetaDirection
    @test direction((u = 0.1, v = -0.2)) isa UVDirection
    @test uv(30°, 0.0) isa UVDirection

    r = region(-10°..10°, 5°; name = :main)
    @test r isa Region
    @test length(r.points) == 5
    @test r.name == :main

    mask = theta_ramp(-10°, -20dB, 10°, -30dB)
    sl = sidelobes(r, mask)
    sb = shaped_beam(r, p -> 1.0; ripple = -0.5dB)
    p = pattern(beam(0°), sb, null(20°), nulls([-30°, 30°]), sl)

    @test length(p.beams) == 1
    @test length(p.shaped_beams) == 1
    @test length(p.null_directions) == 3
    @test length(p.sidelobe_regions) == 1
end

@testset "two-dimensional regions" begin
    shape = Ellipse(0.3, 0.2, (0.0, 0.0))
    beam_region = region(shape; step = 0.1)
    visible = visible_region(shape; step = 0.25, bandpass = 0.05, filtered = true)

    @test beam_region isa Region
    @test visible isa Region
    @test !isempty(beam_region.points)
    @test !isempty(visible.points)
    @test all(p -> p isa UVDirection, visible.points)
    @test all(p -> p.u^2 + p.v^2 < 1, visible.points)

    @test rhombus((0.0, 0.0), 0.2) isa Polygon
    @test triangle((0.0, 0.0); base = 0.2, height = 0.3) isa Polygon
end

@testset "objectives, formulations, excitations, and robustness" begin
    @test ComplexWeights() isa AbstractExcitation
    @test RealWeights() isa AbstractExcitation
    @test ConjugateSymmetricWeights() isa AbstractExcitation
    @test ProgressivePhaseAmplitude() isa AbstractExcitation

    @test LP() isa AbstractFormulation
    @test QP() isa AbstractFormulation
    @test SOCP() isa AbstractFormulation
    @test MILP() isa AbstractFormulation

    r = region(-20°..20°, 10°)
    @test Feasible() isa DirectObjective
    @test MaxAF() isa DirectObjective
    @test MinPower() isa DirectObjective
    @test MinSLL(r) isa DirectObjective
    @test MinIntegratedPower(r) isa DirectObjective
    @test MinL1() isa DirectObjective
    @test MinWeightedL1(ones(4)) isa DirectObjective
    @test MinFieldError(r, 1.0) isa DirectObjective

    @test IterativeReweightedL1() isa SynthesisMethod
    @test IterativeFloorSynthesis(0°) isa SynthesisMethod
    @test MultiPatternReweightedL1() isa SynthesisMethod

    tol = robust(pointing_accuracy = 0.1°, phase_tolerance = 0.5°,
                 amplitude_tolerance = 0.05dB, position_tolerance = (xy = 1e-4, z = 0.0))
    @test tol isa Tolerances
end

@testset "array factor evaluation" begin
    array = uniform_linear_array(3, d = 0.5)
    pts = [ThetaDirection(0.0), ThetaDirection(30°)]
    w = ones(ComplexF64, 3)
    af = array_factor(array, ComplexWeights(), w, pts)
    @test length(af) == length(pts)
    @test af[1] ≈ 3 + 0im

    sym = symmetric_linear_array(3, d = 0.5)
    af_sym = array_factor(sym, ConjugateSymmetricWeights(), ones(ComplexF64, 2), pts)
    @test length(af_sym) == length(pts)
    @test eltype(af_sym) <: Real
end
