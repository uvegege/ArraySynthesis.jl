"""
Abstract supertype for excitation models.
"""
abstract type AbstractExcitation end

"""
    ComplexWeights()

Use independent complex excitation weights for each array element.
"""
struct ComplexWeights <: AbstractExcitation end

"""
    RealWeights()

Use real-valued excitation weights for each array element.
"""
struct RealWeights <: AbstractExcitation end

"""
    ConjugateSymmetricWeights()

Use conjugate-symmetric complex weights on a `SymmetricArray`.
"""
struct ConjugateSymmetricWeights <: AbstractExcitation end

"""
    ProgressivePhaseAmplitude(β)
    ProgressivePhaseAmplitude()

Use real amplitudes with a progressive phase reference `β`.

When `β` is omitted, the reference direction is inferred from the first beam or
shaped-beam region in the pattern.
"""
struct ProgressivePhaseAmplitude <: AbstractExcitation
    β::Union{Real, Tuple{Real, Real}, Vector{Real}, AbstractDirection, Nothing}
end

# Constructor without arguments: β will be extracted from pattern
ProgressivePhaseAmplitude() = ProgressivePhaseAmplitude(nothing)

struct QuantizedAmplitude{L,B} <: AbstractExcitation
    levels::L
    β::B # Like ProgressivePhaseAmplitude
    relative::Bool
end

QuantizedAmplitude(levels::AbstractVector; β = nothing, relative = false) = QuantizedAmplitude(collect(Float64, levels), β, relative)
as_ppa(w::QuantizedAmplitude) = ProgressivePhaseAmplitude(w.β)
