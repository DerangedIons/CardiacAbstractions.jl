"""
    AbstractEPModel

Convenience family root for cardiac-electrophysiology tissue models — inert declarations
of the continuous problem, carrying no mesh, grid, dof handler, or array layout. The real
contract is trait-based: code deciding whether something can be split, or what its
reaction part is, asks [`has_pointwise_reaction_part`](@ref) and the `reaction_*` queries
— never `isa AbstractEPModel` — so a model owned by a third package participates without
re-rooting.
"""
abstract type AbstractEPModel end

"""
    has_pointwise_reaction_part(model) -> Bool

Capability trait: whether `model` has a pointwise reaction part that can be advanced
node-by-node — i.e. whether a [`ReactionDiffusionSplit`](@ref) of it is meaningful.
`false` for arbitrary objects, `true` for the EP family. A foreign model type (an MTK
circuit wrapper, a Dict of subdomain models) opts in by defining this to `true` and
implementing the four `reaction_*` queries: [`reaction_model`](@ref),
[`reaction_solution_symbol`](@ref), [`reaction_state_symbol`](@ref), and
[`reaction_coordinate_system`](@ref).
"""
has_pointwise_reaction_part(model) = false
has_pointwise_reaction_part(::AbstractEPModel) = true

function _reaction_query_error(query::Symbol, model)
    throw(
        ArgumentError(
            "$query is not implemented for $(typeof(model)); a model participating in reaction–diffusion splitting defines `has_pointwise_reaction_part(::$(typeof(model))) = true` plus the four reaction_* queries (reaction_model, reaction_solution_symbol, reaction_state_symbol, reaction_coordinate_system)",
        ),
    )
end

"""
    reaction_model(model)

The pointwise reaction part of `model` — the object satisfying the cell-model contract
that a backend advances at every node. The EP family answers with `model.ion`; foreign
types define their own method.
"""
reaction_model(model::AbstractEPModel) = model.ion
reaction_model(model) = _reaction_query_error(:reaction_model, model)

"""
    reaction_solution_symbol(model) -> Symbol

Published name of the voltage field in solutions of `model` — the symbol solution
accessors resolve. The EP family answers with `model.φ_symbol`.
"""
reaction_solution_symbol(model::AbstractEPModel) = model.φ_symbol
reaction_solution_symbol(model) = _reaction_query_error(:reaction_solution_symbol, model)

"""
    reaction_state_symbol(model) -> Symbol

Published name of the internal-state block of `model` — the non-voltage states as a
whole. The EP family answers with `model.states_symbol`.
"""
reaction_state_symbol(model::AbstractEPModel) = model.states_symbol
reaction_state_symbol(model) = _reaction_query_error(:reaction_state_symbol, model)

"""
    reaction_coordinate_system(model)

The coordinate the reaction part's cell model sees as `x`, or `nothing` for physical
coordinates. The EP family answers with `model.cell_coordinates`; interpreting a
non-`nothing` value (heart axes, transmural depth, …) is the backend's job — this
declaration-level hook is all the base knows about coordinate systems.
"""
reaction_coordinate_system(model::AbstractEPModel) = model.cell_coordinates
reaction_coordinate_system(model) = _reaction_query_error(:reaction_coordinate_system, model)
