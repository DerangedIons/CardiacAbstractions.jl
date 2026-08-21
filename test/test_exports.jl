using CardiacAbstractions
using Test

@testset "export surface is exactly the tree + verbs + traits" begin
    expected = Set(
        [
            :CardiacAbstractions,
            # the type tree
            :AbstractEPModel,
            :MonodomainModel,
            :ReactionDiffusionSplit,
            :AbstractStimulationProtocol,
            :TransmembraneStimulationProtocol,
            :NoStimulationProtocol,
            :AnalyticalTransmembraneStimulationProtocol,
            :AbstractCellModel,
            # the cell-model contract
            :num_states,
            :state_symbols,
            :transmembrane_potential_symbol,
            :transmembrane_potential_index,
            :default_initial_state,
            :cell_rhs!,
            :num_parameters,
            :parameter_names,
            :reaction_rhs!,
            :state_rhs!,
            # capability traits
            :is_cell_model,
            :has_pointwise_reaction_part,
            :reaction_model,
            :reaction_solution_symbol,
            :reaction_state_symbol,
            :reaction_coordinate_system,
            # the verbs
            :semidiscretize,
            :create_initial_condition,
        ]
    )
    @test Set(names(CardiacAbstractions)) == expected
end

@testset "is_active is supported API but deliberately unexported" begin
    @test !(:is_active in names(CardiacAbstractions))
    @test CardiacAbstractions.is_active isa Function
end

@testset "the verbs start as empty generics the base owns" begin
    @test isempty(methods(CardiacAbstractions.semidiscretize))
    @test isempty(methods(CardiacAbstractions.create_initial_condition))
    @test occursin("semidiscretize", string(Base.Docs.doc(CardiacAbstractions.semidiscretize)))
    @test occursin("create_initial_condition", string(Base.Docs.doc(CardiacAbstractions.create_initial_condition)))
end
