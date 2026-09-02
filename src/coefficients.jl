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

struct QuantizedPhase{T, B} <: AbstractExcitation
    amplitude::T
    levels::B
    free_amplitudes::Bool
end

QuantizedPhase(levels::AbstractVector; amplitude = 1.0, free_amplitudes = false) = QuantizedPhase(amplitude, collect(Float64, levels), free_amplitudes)

struct QuantizedWeights <: AbstractExcitation
    constellation::Vector{ComplexF64}
    free_scale::Bool
end

function QuantizedWeights(constellation::AbstractVector{<:Number}; free_scale::Bool = false)
    C = unique(collect(ComplexF64, constellation))
    isempty(C) && error("QuantizedWeights needs at least one constellation point.")
    return QuantizedWeights(C, free_scale)
end

function QuantizedWeights(gains::AbstractVector{<:Real}, phases::AbstractVector{<:Real}; free_scale::Bool = false)
    C = ComplexF64[g * cis(φ) for g in gains for φ in phases]
    return QuantizedWeights(C; free_scale)
end