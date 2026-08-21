"""
    AbstractCellModel

Convenience supertype for cell kinetics — ionic models, mitochondria, signaling, anything
advancing pointwise states. The contract is the *methods*, not this supertype: a type that
cannot subtype it (it already has a supertype in a foreign package) becomes a full citizen
by defining [`is_cell_model`](@ref) to return `true` and implementing the methods below.

# Required

  - `num_states(::Type{M}) -> Int` — number of ODE states (≥ 1); see [`num_states`](@ref).
  - [`default_initial_state`](@ref)`(m) -> AbstractVector` — length `num_states(M)`.
  - [`cell_rhs!`](@ref)`(du, u, x, t, m) -> Nothing` — the kinetics.

# Optional, with derived defaults

  - [`state_symbols`](@ref)`(::Type{M}) -> NTuple{N, Symbol}` — defaults to `(:φₘ, :s1, :s2, …)`.
  - [`transmembrane_potential_symbol`](@ref)`(::Type{M}) -> Symbol` — defaults to `:φₘ`.

[`transmembrane_potential_index`](@ref) is **provided, never implemented**: it is derived
from the two symbol queries, so the published names and the index cannot disagree.

# Optional, no default (probe with `hasmethod`)

  - [`num_parameters`](@ref), [`parameter_names`](@ref) — reflection for atomic models.
  - [`reaction_rhs!`](@ref), [`state_rhs!`](@ref) — sub-split forms.

# Kernel-facing invariants

Every implementation of `cell_rhs!` writes every slot of `du` before any slot is
read-modified (never `+=` into an unwritten slot) and computes in `eltype(u)` — wrap
numeric literals as `T(0.5)` so a `Float32` run stays `Float32` end to end. Type-level
queries return constants so they fold at compile time and no `Symbol` survives into a hot
kernel.
"""
abstract type AbstractCellModel end

"""
    is_cell_model(m) -> Bool

Whether `m` satisfies the cell-model contract. `true` for any [`AbstractCellModel`](@ref);
a foreign type opts in — without re-rooting under the base's supertype — by defining a
method returning `true` and implementing the required contract methods.
"""
is_cell_model(m) = m isa AbstractCellModel

"""
    num_states(::Type{M}) -> Int
    num_states(m) -> Int

Number of ODE states of the cell model. Required, `≥ 1`, and a property of the model
*type* so it folds to a compile-time constant; the instance form forwards to the type
form.
"""
num_states(m) = num_states(typeof(m))

# The `::Type` fallback turns a missing implementation into an actionable error. Without
# it, calling `num_states` on an unimplemented model *type* would match the instance
# forwarder and recurse through `typeof(M) === DataType` into a silent StackOverflowError.
function num_states(::Type{T}) where {T}
    throw(
        ArgumentError(
            "$T does not implement the cell-model contract: define `CardiacAbstractions.num_states(::Type{$T}) = n` along with `default_initial_state` and `cell_rhs!` (see the AbstractCellModel docstring)",
        ),
    )
end

"""
    state_symbols(::Type{M}) -> NTuple{N, Symbol}
    state_symbols(m)

Published names of the state slots, in storage order. Defaults to `(:φₘ, :s1, :s2, …)`;
override with a literal tuple to name states. The voltage slot is identified by
[`transmembrane_potential_symbol`](@ref) and its index is always derived — never stated —
via [`transmembrane_potential_index`](@ref).
"""
state_symbols(m) = state_symbols(typeof(m))

@inline function state_symbols(::Type{T}) where {T}
    return ntuple(i -> i == 1 ? :φₘ : Symbol(:s, i - 1), Val(num_states(T)))
end

"""
    transmembrane_potential_symbol(::Type{M}) -> Symbol
    transmembrane_potential_symbol(m)

The *role* marker naming which entry of [`state_symbols`](@ref) is the transmembrane
potential. Defaults to `:φₘ`. A role, not a user-facing name: a model whose voltage state
is called `:V` overrides this to `:V` and keeps its own naming.
"""
transmembrane_potential_symbol(m) = transmembrane_potential_symbol(typeof(m))
transmembrane_potential_symbol(::Type{T}) where {T} = :φₘ

# Recursive tuple search instead of `findfirst`: every recursion level peels one element
# off the tuple *type*, so inference unrolls the search structurally and the result folds
# to a compile-time constant by inlining alone — no `@generated` (deliberately, for the
# world-age reason documented in Thunderbolt).
@inline _symbol_index(::Tuple{}, target::Symbol, i::Int) = nothing
@inline function _symbol_index(syms::Tuple, target::Symbol, i::Int)
    return first(syms) === target ? i : _symbol_index(Base.tail(syms), target, i + 1)
end

"""
    transmembrane_potential_index(::Type{M}) -> Int
    transmembrane_potential_index(m)

Index of the transmembrane potential in the state vector. Provided and **derived** from
[`state_symbols`](@ref) and [`transmembrane_potential_symbol`](@ref) — do not implement
it, so names and index cannot disagree. Folds to a compile-time constant.
"""
@inline transmembrane_potential_index(m) = transmembrane_potential_index(typeof(m))

@inline function transmembrane_potential_index(::Type{T}) where {T}
    sym = transmembrane_potential_symbol(T)
    idx = _symbol_index(state_symbols(T), sym, 1)
    idx === nothing && throw(
        ArgumentError(
            "transmembrane_potential_symbol(::Type{$T}) = :$sym does not appear in state_symbols(::Type{$T}) = $(state_symbols(T))",
        ),
    )
    return idx
end

"""
    default_initial_state(m) -> AbstractVector

Initial condition of length `num_states(typeof(m))`, in the model's own units. Required.
Plain CPU memory; backends relocate it to wherever their semidiscrete function lives (see
[`create_initial_condition`](@ref)).
"""
function default_initial_state end

"""
    cell_rhs!(du, u, x, t, m) -> Nothing

The cell kinetics: write d(state)/dt into `du` given states `u`, position `x`, and time
`t`. Required.

  - `x` is the coordinate the surrounding tissue model publishes — physical position
    unless the model declares `cell_coordinates` — and may be `nothing` for
    coordinate-free models, which must not dereference it.
  - **Write every slot of `du`** before any slot is read-modified: `du` arrives
    uninitialized, so `du[i] += v` as a slot's first write reads garbage.
  - **Compute in `eltype(u)`** — wrap numeric literals in `T(…)` so a `Float32` run stays
    `Float32`.
  - Spatial parameter heterogeneity never widens this signature: backends curry the tissue
    model's `overrides` into the cell model at [`semidiscretize`](@ref) time.
  - The tissue model owns the stimulus. A cell model used on the tissue path should carry
    none of its own; backends may warn when one does.
"""
function cell_rhs!(du, u, x, t, m)
    # Catch-all so a missing (or mis-signatured) implementation fails with an actionable
    # message instead of a MethodError. `cell_rhs!` is required, so `hasmethod` is never
    # the probe for it — unlike the optional methods below, which stay bare on purpose.
    throw(
        ArgumentError(
            "$(typeof(m)) does not implement `CardiacAbstractions.cell_rhs!(du, u, x, t, m)`, required by the cell-model contract (see the AbstractCellModel docstring)",
        ),
    )
end

"""
    num_parameters(m) -> Int

Optional reflection: number of parameters of an *atomic* model, for fitting and parameter
studies. A composite model leaves it undefined so `hasmethod` correctly reports the
missing capability — no fallback is provided on purpose.
"""
function num_parameters end

"""
    parameter_names(m) -> NTuple{N, Symbol}

Optional reflection: parameter names of an atomic model, in the order a parameter vector
uses. Undefined for composites; probe with `hasmethod`.
"""
function parameter_names end

"""
    reaction_rhs!(dφ, φ, s, x, t, m) -> Nothing

Optional sub-split form: the voltage equation alone, for backends splitting within the
cell model. Same invariants as [`cell_rhs!`](@ref). Probe with `hasmethod`.
"""
function reaction_rhs! end

"""
    state_rhs!(ds, φ, s, x, t, m) -> Nothing

Optional sub-split form: the non-voltage state equations alone, for backends splitting
within the cell model. Same invariants as [`cell_rhs!`](@ref). Probe with `hasmethod`.
"""
function state_rhs! end

# DiffEq-style convenience on the coordinate-free path only. Any non-`nothing` `p` keeps
# its semantics with the model's owner (spatial context, curried overrides, …), so no
# method is provided for it here.
@inline (m::AbstractCellModel)(du, u, ::Nothing, t) = cell_rhs!(du, u, nothing, t, m)
