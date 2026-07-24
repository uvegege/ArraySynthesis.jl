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
        "Home" => "index.md",
        "User guide" => [
            "Overview" => "guide/workflow.md",
            "Pattern and regions" => "guide/patterns.md",
            "Objectives" => "guide/objectives.md",
            "Excitations" => "guide/excitations.md",
            "Formulations" => "guide/formulations.md",
            "Robust synthesis" => "guide/robustness.md",
        ],
        "Guided examples" => [
            "Flat-top linear beam"          => "examples/flat_top_linear.md",
            "Ramp sidelobe mask"            => "examples/narrow_beam_ramp_sll.md",
            "Two-level sidelobe mask"       => "examples/nonuniform_sll_steps.md",
            "Planar shaped beam"            => "examples/planar_shaped_beam.md",
            "Power and null constraints"    => "examples/power_nulls_compare.md",
            "Sparse multipattern thinning"  => "examples/sparse_multipattern.md",
        ],
        "Theory" => "theory.md",
        "Implementation notes" => "implementation.md",
        "API Reference" => "api.md",
    ],
    checkdocs = :none,
    remotes   = nothing,
    warnonly  = [:docs_block, :cross_references],
)

deploydocs(
    repo = "github.com/uvegege/ArraySynthesis.jl.git",
    devbranch = "main",
    forcepush = true
)