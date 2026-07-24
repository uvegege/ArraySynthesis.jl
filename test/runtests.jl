using ArraySynthesis
using ArraySynthesis: °, dB
using JuMP
using Test

const AS = ArraySynthesis

constraint_count(model) = num_constraints(model; count_variable_in_set_constraints = true)

@testset "array geometry" begin
    ula = uniform_linear_array(4, d = 0.5)
    @test ula isa ArrayGeometry
    @test ula.dim == 1
    @test size(ula.positions) == (3, 4)
    @test ula.positions[1, :] ≈ [-0.75, -0.25, 0.25, 0.75]
    @test all(iszero, ula.positions[2:3, :])

    yline = linear_array([-1.0, 0.0, 1.0], axis = :y)
    @test yline.positions[2, :] == [-1.0, 0.0, 1.0]
    @test all(iszero, yline.positions[[1, 3], :])

    planar = planar_array(2, 3, dx = 0.5, dy = 0.25)
    @test planar isa ArrayGeometry
    @test planar.dim == 2
    @test size(planar.positions) == (3, 6)
    @test length(unique(planar.positions[1, :])) == 2
    @test length(unique(planar.positions[2, :])) == 3

    @test size(triangular_array(2, 2).positions, 2) == 4
    @test size(hexagonal_array(1).positions, 2) == 7
    @test size(circular_array([0.0, 1.0], [1, 4]).positions, 2) == 5
end

@testset "symmetric arrays" begin
    sym = symmetric_linear_array(5, d = 0.5)
    @test sym isa SymmetricArray
    @test sym.dim == 1
    @test size(sym.positions, 2) == 3
    @test count(p -> is_origin(p), eachcol(sym.positions)) == 1

    full = materialize(sym)
    @test full isa ArrayGeometry
    @test size(full.positions, 2) == 5
    @test sort(full.positions[1, :]) ≈ [-1.0, -0.5, 0.0, 0.5, 1.0]

    reps, dim = symmetrize(full.positions)
    @test dim == 1
    @test size(reps, 2) == size(sym.positions, 2)
    @test all(col -> col[1] >= 0, eachcol(reps))
end

@testset "directions and regions" begin
    @test direction(10°) == ThetaDirection(10°)
    @test direction((u = 0.2, v = -0.1)) == UVDirection(0.2, -0.1)
    @test uv(30°, 0.0).u ≈ sin(30°)
    @test uv(30°, 0.0).v ≈ 0.0

    r = region(-10°..10°, 5°; name = :main)
    @test r.name == :main
    @test [p.θ for p in r.points] ≈ collect(-10°:5°:10°)

    pieces = outside([-20°..(-10°), 10°..20°])
    @test length(pieces) == 3
    @test pieces[1].a ≈ -90°
    @test pieces[1].b ≈ -20°
    @test pieces[2].a ≈ -10°
    @test pieces[2].b ≈ 10°
    @test pieces[3].a ≈ 20°
    @test pieces[3].b ≈ 90°

    joined = reduce(join_regions, region.(pieces, 10°))
    @test length(joined.points) == sum(length, (region(p, 10°).points for p in pieces))
end

@testset "two-dimensional shapes" begin
    circle = Circle(0.3, (0.0, 0.0))
    ellipse = Ellipse(0.3, 0.2, (0.0, 0.0))
    rect = Rectangle(0.2, 0.1, (0.0, 0.0))

    @test AS.isinside(circle, 0.0, 0.0)
    @test !AS.isinside(circle, 0.4, 0.0)
    @test AS.isinside(ellipse, 0.0, 0.0)
    @test AS.isinside(rect, 0.1, 0.05)

    beam_region = region(ellipse; step = 0.1)
    visible = visible_region(ellipse; step = 0.25, bandpass = 0.05, filtered = true)
    @test !isempty(beam_region.points)
    @test !isempty(visible.points)
    @test all(p -> p isa UVDirection, visible.points)
    @test all(p -> p.u^2 + p.v^2 < 1, visible.points)
    @test all(p -> !AS.isinside(ellipse, p.u, p.v), visible.points)

    @test rhombus((0.0, 0.0), 0.2) isa Polygon
    @test triangle((0.0, 0.0); base = 0.2, height = 0.3) isa Polygon
end

@testset "pattern assembly" begin
    main = region(-10°..10°, 5°)
    sl_region = region(20°..40°, 10°)
    mask = theta_ramp(20°, -20dB, 40°, -30dB)

    p = pattern(
        sidelobes(sl_region, mask),
        nulls([-40°, 40°]),
        shaped_beam(main, x -> 1.0; ripple = -0.5dB),
        beam(0°),
        null(50°),
    )

    @test length(p.beams) == 1
    @test p.beams[1].gain == 1.0
    @test length(p.shaped_beams) == 1
    @test p.shaped_beams[1].region === main
    @test length(p.null_directions) == 3
    @test length(p.sidelobe_regions) == 1
    @test p.sidelobe_regions[1].upper[1] ≈ -20dB
    @test p.sidelobe_regions[1].upper[end] ≈ -30dB
end

@testset "array factor evaluation" begin
    array = uniform_linear_array(3, d = 0.5)
    pts = [ThetaDirection(0.0), ThetaDirection(30°)]
    w = ones(ComplexF64, 3)

    af = array_factor(array, ComplexWeights(), w, pts)
    @test length(af) == 2
    @test af[1] ≈ 3 + 0im

    real_af = array_factor(array, RealWeights(), ones(3), pts)
    @test real_af ≈ af

    sym = symmetric_linear_array(3, d = 0.5)
    sym_af = array_factor(sym, ConjugateSymmetricWeights(), ones(ComplexF64, 2), pts)
    @test length(sym_af) == length(pts)
    @test eltype(sym_af) <: Real
end

@testset "model variables and compatibility errors" begin
    model = Model()
    array = uniform_linear_array(4, d = 0.5)
    sym = symmetric_linear_array(4, d = 0.5)

    complex_vars = AS.variables!(model, array, ComplexWeights(), LP())
    @test length(complex_vars.w_re) == 4
    @test length(complex_vars.w_im) == 4

    real_vars = AS.variables!(Model(), array, RealWeights(), LP())
    @test length(real_vars.w_re) == 4
    @test isempty(real_vars.w_im)

    amp_vars = AS.variables!(Model(), array, ProgressivePhaseAmplitude(0°), LP())
    @test length(amp_vars.a) == 4

    @test_throws ErrorException AS.variables!(Model(), array, ConjugateSymmetricWeights(), LP())
    @test_throws ErrorException AS.variables!(Model(), sym, ComplexWeights(), LP())
    @test_throws ErrorException AS.variables!(Model(), sym, RealWeights(), LP())

    @test_throws ErrorException AS.check_formulation(MinPower(), LP())
    @test_throws ErrorException AS.check_formulation(MinIntegratedPower(region(-10°..10°, 5°)), LP())
    @test_throws ErrorException AS.check_formulation(MinFieldError(region(-10°..10°, 5°), 1.0), LP())
    @test try
        AS.check_formulation(MinPower(), QP())
        true
    catch
        false
    end
end

@testset "objective and constraint construction" begin
    array = symmetric_linear_array(4, d = 0.5)
    weights = ConjugateSymmetricWeights()
    formulation = LP()
    r = region(-30°..30°, 30°)
    p = pattern(beam(0°), sidelobes(r, -20dB))

    model = Model()
    vars = AS.variables!(model, array, weights, formulation)
    AS.constraints!(model, p, array, weights, vars, formulation)
    @test num_variables(model) == 4
    @test constraint_count(model) > 0

    model = Model()
    vars = AS.variables!(model, array, weights, formulation)
    aux = AS.objective!(model, MaxAF(), p, array, weights, vars, formulation)
    @test aux === nothing
    @test objective_sense(model) == MOI.MAX_SENSE

    model = Model()
    vars = AS.variables!(model, array, weights, formulation)
    t = AS.objective!(model, MinSLL(r), p, array, weights, vars, formulation)
    @test length(t) == 1
    @test objective_sense(model) == MOI.MIN_SENSE

    model = Model()
    vars = AS.variables!(model, array, weights, SOCP())
    AS.constraints!(model, pattern(shaped_beam(r, 1.0, ripple = -0.5dB)), array, weights, vars, SOCP(), nothing)
    @test constraint_count(model) > 0
end
