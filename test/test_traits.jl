using CardiacAbstractions
using Test

struct StubIon <: AbstractCellModel end
CardiacAbstractions.num_states(::Type{StubIon}) = 2
CardiacAbstractions.default_initial_state(::StubIon) = [0.0, 0.0]
function CardiacAbstractions.cell_rhs!(du, u, x, t, m::StubIon)
    du[1] = -u[1]
    du[2] = -u[2]
    return nothing
end

# ---------------------------------------------------------------------------
# A third package's model and ion: they subtype NEITHER abstract root and opt in
# purely through the trait grammar and the contract methods (design principle 2).
# ---------------------------------------------------------------------------
struct ForeignIon end
CardiacAbstractions.is_cell_model(::ForeignIon) = true
CardiacAbstractions.num_states(::Type{ForeignIon}) = 2
CardiacAbstractions.state_symbols(::Type{ForeignIon}) = (:φₘ, :w)
CardiacAbstractions.default_initial_state(::ForeignIon) = [0.0, 0.0]
function CardiacAbstractions.cell_rhs!(du, u, x, t, m::ForeignIon)
    du[1] = -u[1]
    du[2] = u[1] - u[2]
    return nothing
end

struct CircuitWrapper{C}
    cell::C
end
CardiacAbstractions.has_pointwise_reaction_part(::CircuitWrapper) = true
CardiacAbstractions.reaction_model(w::CircuitWrapper) = w.cell
CardiacAbstractions.reaction_solution_symbol(::CircuitWrapper) = :φₘ
CardiacAbstractions.reaction_state_symbol(::CircuitWrapper) = :states
CardiacAbstractions.reaction_coordinate_system(::CircuitWrapper) = nothing

idxval(m) = Val(transmembrane_potential_index(m))

@testset "the EP family answers the traits through its fields" begin
    model = MonodomainModel(; κ = 1.0e-3, ion = StubIon())
    @test has_pointwise_reaction_part(model)
    @test reaction_model(model) === model.ion
    @test reaction_solution_symbol(model) === :φₘ
    @test reaction_state_symbol(model) === :states
    @test reaction_coordinate_system(model) === nothing

    renamed = MonodomainModel(;
        κ = 1.0e-3,
        ion = StubIon(),
        φ_symbol = :V,
        states_symbol = :gates,
        cell_coordinates = :transmural,
    )
    @test reaction_solution_symbol(renamed) === :V
    @test reaction_state_symbol(renamed) === :gates
    @test reaction_coordinate_system(renamed) === :transmural
end

@testset "arbitrary objects are not splittable" begin
    @test !has_pointwise_reaction_part(42)
    @test !has_pointwise_reaction_part("tissue")
    @test !has_pointwise_reaction_part(nothing)
end

@testset "a foreign type is a full citizen without re-rooting" begin
    @test !(ForeignIon <: AbstractCellModel)
    @test !(CircuitWrapper <: AbstractEPModel)
    ion = ForeignIon()
    wrapper = CircuitWrapper(ion)
    @test is_cell_model(ion)
    @test has_pointwise_reaction_part(wrapper)
    @test reaction_model(wrapper) === ion
    @test reaction_solution_symbol(wrapper) === :φₘ
    @test reaction_state_symbol(wrapper) === :states
    @test reaction_coordinate_system(wrapper) === nothing
    @test transmembrane_potential_index(ion) == 1
    @test @inferred(idxval(ion)) === Val(1)
    split = ReactionDiffusionSplit(wrapper)    # M is unconstrained on purpose
    @test split.model === wrapper
end

@testset "reaction_* queries teach the opt-in when unimplemented" begin
    @test_throws ArgumentError reaction_model(42)
    @test_throws "has_pointwise_reaction_part" reaction_model(42)
    @test_throws ArgumentError reaction_solution_symbol(42)
    @test_throws ArgumentError reaction_state_symbol(42)
    @test_throws ArgumentError reaction_coordinate_system(42)
end
