using CardiacAbstractions
using Test
using SafeTestsets
using JET

@testset "CardiacAbstractions.jl" begin
    @testset "Code linting (JET.jl)" begin
        # `target_modules` rather than `target_defined_modules`: JET 0.12 dropped the
        # latter configuration name. Same restriction — report only problems inside
        # CardiacAbstractions' own module context.
        JET.test_package(CardiacAbstractions; target_modules = (CardiacAbstractions,))
    end

    # `test_exports.jl` asserts the verbs' method tables are still empty, so it must run
    # before `test_verbs.jl` defines the fake backend's methods (method tables are
    # global; SafeTestsets does not undo definitions).
    @safetestset "Export surface" begin
        include("test_exports.jl")
    end
    @safetestset "Cell-model contract" begin
        include("test_cell_interface.jl")
    end
    @safetestset "Stimulation vocabulary" begin
        include("test_stimulation.jl")
    end
    @safetestset "Capability traits" begin
        include("test_traits.jl")
    end
    @safetestset "MonodomainModel" begin
        include("test_monodomain.jl")
    end
    @safetestset "ReactionDiffusionSplit" begin
        include("test_split.jl")
    end
    @safetestset "Verbs and canonical programs" begin
        include("test_verbs.jl")
    end
end
