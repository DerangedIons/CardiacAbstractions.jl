using CardiacAbstractions
using Test

struct SplitIon <: AbstractCellModel end
CardiacAbstractions.num_states(::Type{SplitIon}) = 2

@testset "a pure annotation" begin
    model = MonodomainModel(; κ = 1.0e-3, ion = SplitIon())
    split = ReactionDiffusionSplit(model)
    @test split.model === model
    # overrides live on the model now, not on the split
    @test fieldnames(ReactionDiffusionSplit) == (:model,)
    # no abstract split hierarchy until a second split annotation exists
    @test supertype(ReactionDiffusionSplit) === Any
end

@testset "the model slot is deliberately unconstrained" begin
    # trait-checked at semidiscretize, not type-constrained here
    @test ReactionDiffusionSplit(42).model === 42
end

@testset "wrapping is free" begin
    @test isbitstype(typeof(ReactionDiffusionSplit(NoStimulationProtocol())))
    model = MonodomainModel(; κ = 1.0e-3, ion = SplitIon())
    @test @inferred(ReactionDiffusionSplit(model)) isa ReactionDiffusionSplit
end
