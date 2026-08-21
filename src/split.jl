"""
    ReactionDiffusionSplit(model)

Annotate `model` for a reaction–diffusion operator split: the pointwise reaction part is
advanced node-by-node, the diffusion operator globally. A pure annotation saying *how* to
split and nothing else — no work happens until [`semidiscretize`](@ref), and
`semidiscretize(model, disc, geometry)` without the wrapper is the monolithic path (one
coupled right-hand side).

The model slot is deliberately unconstrained: `semidiscretize` trait-checks
[`has_pointwise_reaction_part`](@ref) instead of a supertype, so foreign model types can
be split. There is no abstract split hierarchy — none exists until a second split
annotation does.
"""
struct ReactionDiffusionSplit{M}
    model::M
end
