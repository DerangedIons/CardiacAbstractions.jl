"""
    semidiscretize(model_or_split, discretization, geometry)

Turn a continuous declaration into a spatially-discrete, time-continuous function ready
for a time integrator. Owned by CardiacAbstractions as an empty generic function,
CommonSolve-style: the base fixes the argument grammar and the promise, and each backend
adds methods for its own discretization and geometry types — two backends defining
methods of the *same* function cannot drift on its meaning.

# Contract

  - Passing a model is the monolithic path (one coupled right-hand side); wrapping it in
    [`ReactionDiffusionSplit`](@ref) is the declarative request for an operator-splitting
    function whose reaction/diffusion structure the downstream integrator consumes.
  - Backends decide splittability with [`has_pointwise_reaction_part`](@ref) and read the
    reaction part through the `reaction_*` queries — never `isa`.
  - Backends must reject coefficient or override values they cannot interpret *here*, at
    semidiscretization time — never at solve time. The portable subset every backend
    honors is `Number` (homogeneous) and callables (analytic fields).

See also: [`create_initial_condition`](@ref).
"""
function semidiscretize end

"""
    create_initial_condition(f) -> u₀

Allocate the solution vector for a semidiscrete function produced by
[`semidiscretize`](@ref), filled from the reaction part's
[`default_initial_state`](@ref), on the device and in the layout `f` lives on. Owned by
CardiacAbstractions as an empty generic function; backends add methods for their own
semidiscrete function types.
"""
function create_initial_condition end
