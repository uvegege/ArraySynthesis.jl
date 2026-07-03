abstract type AbstractExcitation end

struct ComplexWeights <: AbstractExcitation end
struct RealWeights <: AbstractExcitation end
struct ConjugateSymmetricWeights <: AbstractExcitation end

struct ProgressivePhaseAmplitude <: AbstractExcitation
    β::Union{Real, Tuple{Real, Real}, Vector{Real}, AbstractDirection, Nothing}
end

# Constructor without arguments: β will be extracted from pattern
ProgressivePhaseAmplitude() = ProgressivePhaseAmplitude(nothing)

