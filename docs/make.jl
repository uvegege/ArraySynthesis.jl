using Documenter
using Literate
using ArraySynthesis

literate_src = joinpath(@__DIR__, "literate")
literate_dst = joinpath(@__DIR__, "src", "examples")
mkpath(literate_dst)

for file in sort(filter(endswith(".jl"), readdir(literate_src; join = true)))
    Literate.markdown(file, literate_dst; documenter = false, execute = false)
end

makedocs(
    sitename = "ArraySynthesis.jl",
    modules  = [ArraySynthesis],
    format   = Documenter.HTML(),
    pages    = [
        "Home"    => "index.md",
        "Theory"  => "theory.md",
        "Objectives" => [
            "Overview"                     => "guide/objectives.md",
            "Feasible"                     => "guide/objectives/feasible.md",
            "MinSLL"                       => "guide/objectives/minsll.md",
            "MaxGain"                      => "guide/objectives/maxgain.md",
            "MinL1"                        => "guide/objectives/minl1.md",
            "MinWeightedL1"                => "guide/objectives/minweightedl1.md",
            "IterativeReweightedL1"        => "guide/objectives/iterativereweightedl1.md",
            "MultiPatternReweightedL1"     => "guide/objectives/multipatternreweightedl1.md",
            "MinPower"                     => "guide/objectives/minpower.md",
            "MinIntegratedPower"           => "guide/objectives/minintegratedpower.md",
            "MinFieldError"                => "guide/objectives/minfielderror.md",
            "IterativePatternLeastSquares" => "guide/objectives/iterativepatternls.md",
            "IterativeFloorSynthesis"      => "guide/objectives/iterativefloor.md",
            "MaxDirectivity"               => "guide/objectives/maxdirectivity.md",
        ],
        "Excitations" => [
            "Overview"                  => "guide/excitations.md",
            "ComplexWeights"            => "guide/excitations/complexweights.md",
            "RealWeights"               => "guide/excitations/realweights.md",
            "ConjugateSymmetricWeights" => "guide/excitations/conjugatesymmetric.md",
            "ProgressivePhaseAmplitude" => "guide/excitations/progressivephase.md",
        ],
        "Formulations" => [
            "Overview" => "guide/formulations.md",
            "LP"   => "guide/formulations/lp.md",
            "QP"   => "guide/formulations/qp.md",
            "SOCP" => "guide/formulations/socp.md",
        ],
        "Robustness" => "guide/robustness.md",
        "Guided examples" => [
            "Flat-top linear beam"          => "examples/flat_top_linear.md",
            "Ramp sidelobe mask"            => "examples/narrow_beam_ramp_sll.md",
            "Two-level sidelobe mask"       => "examples/nonuniform_sll_steps.md",
            "Planar shaped beam"            => "examples/planar_shaped_beam.md",
            "Power and null constraints"    => "examples/power_nulls_compare.md",
            "Sparse multipattern thinning"  => "examples/sparse_multipattern.md",
        ],
        "Implementation notes" => "implementation.md",
        "API"     => "api.md",
    ],
    checkdocs = :none,
    remotes   = nothing,
    warnonly  = [:docs_block, :cross_references],
)
