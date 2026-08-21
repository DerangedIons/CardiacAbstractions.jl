"""
    AbstractStimulationProtocol

Root of the stimulation vocabulary: an externally applied stimulation current, declared
at the continuum level as a callable `(x, t) -> Iₛₜᵢₘ`. Deliberately wider than the transmembrane
family — split intra-/extracellular stimulation re-enters through this seam when bidomain
declarations join the base.

# Sign convention

Stimulus values are in the PDE-form convention: **positive is depolarizing**. Stated once,
here, and inherited by every backend.
"""
abstract type AbstractStimulationProtocol end

"""
    TransmembraneStimulationProtocol <: AbstractStimulationProtocol

The physical statement `Iₛₜᵢₘ,ᵢ = Iₛₜᵢₘ,ₑ` — a stimulus applied equally to the
intracellular and extracellular spaces, which is what enters a monodomain right-hand side
as a transmembrane source term. Abstract: the concrete analytic shape is
[`AnalyticalTransmembraneStimulationProtocol`](@ref), and the absence of stimulation is
[`NoStimulationProtocol`](@ref).
"""
abstract type TransmembraneStimulationProtocol <: AbstractStimulationProtocol end

"""
    NoStimulationProtocol()

The absence of an applied stimulus. Fieldless on purpose: dispatching on it deletes the
stimulus term from a right-hand side entirely — a strong zero — rather than evaluating a
function that returns zero. Its functor returns `false`, the `Bool` zero that promotes to
whatever it lands in.
"""
struct NoStimulationProtocol <: TransmembraneStimulationProtocol end

@inline (::NoStimulationProtocol)(x, t) = false

_normalize_intervals(::Nothing) = nothing

function _normalize_intervals(intervals)
    normalized = Tuple(
        map(intervals) do iv
            length(iv) == 2 || throw(
                ArgumentError("each nonzero interval must be a (t₀, t₁) pair, got $iv")
            )
            t₀, t₁ = iv
            t₀ <= t₁ || throw(
                ArgumentError("nonzero interval $iv is empty: t₀ must not exceed t₁")
            )
            return (t₀, t₁)
        end,
    )
    isempty(normalized) && throw(
        ArgumentError(
            "nonzero_intervals is empty, which would silence the protocol at every time; pass `nothing` to evaluate it unconditionally",
        ),
    )
    return normalized
end

"""
    AnalyticalTransmembraneStimulationProtocol(; f, nonzero_intervals = nothing)

A transmembrane stimulus given by a plain callable `f(x, t) -> Iₛₜᵢₘ`, positive =
depolarizing. The declaration is continuum-level: a structured-grid backend evaluates `f`
pointwise, an FEM backend lowers it into its coefficient machinery at
[`semidiscretize`](@ref) time.

`nonzero_intervals` is an optional collection of `(t₀, t₁)` windows outside which `f` is
known to vanish — a sparsity-in-time hint every backend may exploit through
[`CardiacAbstractions.is_active`](@ref), worth having when a 2 ms pulse rides on a 500 ms
simulation. It is normalized to a `Tuple` of pairs, never a `Vector`, so the protocol is
`isbits` exactly when `f` is: a plain function or a small immutable functor rides into a
GPU kernel by value, a closure over a boxed variable does not.

### Examples

```julia
# 2 ms depolarizing pulse over the left 1.5 length units of the domain
stim = AnalyticalTransmembraneStimulationProtocol(;
    f = (x, t) -> (t <= 2.0 && x[1] <= 1.5) ? 50.0 : 0.0,
    nonzero_intervals = ((0.0, 2.0),),
)
```
"""
struct AnalyticalTransmembraneStimulationProtocol{F, W} <: TransmembraneStimulationProtocol
    f::F
    nonzero_intervals::W

    # Keyword-only inner constructor: defining it suppresses the generated positional
    # constructors, so no construction path bypasses interval validation.
    function AnalyticalTransmembraneStimulationProtocol(; f, nonzero_intervals = nothing)
        intervals = _normalize_intervals(nonzero_intervals)
        return new{typeof(f), typeof(intervals)}(f, intervals)
    end
end

@inline (p::AnalyticalTransmembraneStimulationProtocol)(x, t) = p.f(x, t)

"""
    is_active(protocol, t) -> Bool

Whether `protocol` can be nonzero at time `t` — the consumer of `nonzero_intervals`.
Conservative: a protocol without declared windows, or a foreign protocol type, is always
reported active. Backends use it to skip stimulus evaluation, never to decide the value.

Unexported but supported API — call it as `CardiacAbstractions.is_active`.
"""
is_active(::AbstractStimulationProtocol, t) = true
is_active(::NoStimulationProtocol, t) = false

# The no-windows case is a type-parameter dispatch, not a runtime `=== nothing` branch.
is_active(::AnalyticalTransmembraneStimulationProtocol{F, Nothing}, t) where {F} = true

function is_active(p::AnalyticalTransmembraneStimulationProtocol, t)
    for (t₀, t₁) in p.nonzero_intervals
        t₀ <= t <= t₁ && return true
    end
    return false
end
