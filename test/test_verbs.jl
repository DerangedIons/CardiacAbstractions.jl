using CardiacAbstractions
using Test

# ---------------------------------------------------------------------------
# A complete tiny backend, defined in-test: enough to prove the verb grammar
# composes and that a backend drives everything through the traits — never `isa`.
# ---------------------------------------------------------------------------
struct FakeDiscretization end

struct FakeGrid{B, N}
    bounds::B
    dims::NTuple{N, Int}
end

struct FakeSemidiscreteFunction{S, G}
    split::S
    geo::G
    n::Int
end

function CardiacAbstractions.semidiscretize(
        split::ReactionDiffusionSplit, ::FakeDiscretization, geo::FakeGrid
    )
    m = split.model
    has_pointwise_reaction_part(m) || throw(
        ArgumentError("$(typeof(m)) has no pointwise reaction part; it cannot be split")
    )
    ion = reaction_model(m)
    n = prod(geo.dims) * num_states(ion)
    return FakeSemidiscreteFunction(split, geo, n)
end

function CardiacAbstractions.semidiscretize(
        model::MonodomainModel, disc::FakeDiscretization, geo::FakeGrid
    )
    # the monolithic spelling of the same grammar
    return CardiacAbstractions.semidiscretize(ReactionDiffusionSplit(model), disc, geo)
end

function CardiacAbstractions.create_initial_condition(f::FakeSemidiscreteFunction)
    u0 = default_initial_state(reaction_model(f.split.model))
    return repeat(u0, prod(f.geo.dims))
end

# the two backend lines of §3 program 1 run byte-identical under these aliases
const CartesianGrid = FakeGrid
const FiniteDifferenceDiscretization = FakeDiscretization

struct VerbIon <: AbstractCellModel end
CardiacAbstractions.num_states(::Type{VerbIon}) = 2
CardiacAbstractions.default_initial_state(::VerbIon) = [0.0, 1.0]
function CardiacAbstractions.cell_rhs!(du, u, x, t, m::VerbIon)
    du[1] = -u[1]
    du[2] = -u[2]
    return nothing
end

struct ForeignIon end
CardiacAbstractions.is_cell_model(::ForeignIon) = true
CardiacAbstractions.num_states(::Type{ForeignIon}) = 3
CardiacAbstractions.default_initial_state(::ForeignIon) = [0.0, 0.0, 0.0]

struct CircuitWrapper{C}
    cell::C
end
CardiacAbstractions.has_pointwise_reaction_part(::CircuitWrapper) = true
CardiacAbstractions.reaction_model(w::CircuitWrapper) = w.cell
CardiacAbstractions.reaction_solution_symbol(::CircuitWrapper) = :φₘ
CardiacAbstractions.reaction_state_symbol(::CircuitWrapper) = :states
CardiacAbstractions.reaction_coordinate_system(::CircuitWrapper) = nothing

@testset "§3 program 1 composes end-to-end at the declaration level" begin
    model = MonodomainModel(;
        κ = (0.133, 0.0176, 0.0176),
        χ = 140.0,
        Cₘ = 0.01,
        ion = VerbIon(),
        stim = AnalyticalTransmembraneStimulationProtocol(;
            f = (x, t) -> (t ≤ 2.0 && x[1] ≤ 1.5) ? 50.0 : 0.0,
            nonzero_intervals = ((0.0, 2.0),),
        ),
    )
    split = ReactionDiffusionSplit(model)

    geo = CartesianGrid(((0.0, 20.0), (0.0, 7.0), (0.0, 3.0)), (100, 35, 15))  # ← backend line 1
    f = semidiscretize(split, FiniteDifferenceDiscretization(), geo)           # ← backend line 2

    @test f isa FakeSemidiscreteFunction
    @test f.n == 100 * 35 * 15 * 2

    u₀ = create_initial_condition(f)
    @test length(u₀) == f.n     # the verb's promise: sized for f, allocated where f lives
    @test u₀[1:2] == [0.0, 1.0]
end

@testset "monolithic and split spellings both dispatch" begin
    model = MonodomainModel(; κ = 1.0e-3, ion = VerbIon())
    geo = FakeGrid(((0.0, 1.0),), (8,))
    f_mono = semidiscretize(model, FakeDiscretization(), geo)
    f_split = semidiscretize(ReactionDiffusionSplit(model), FakeDiscretization(), geo)
    @test f_mono.n == f_split.n == 8 * 2
end

@testset "a foreign model splits through the traits alone" begin
    wrapper = CircuitWrapper(ForeignIon())
    geo = FakeGrid(((0.0, 1.0),), (4,))
    f = semidiscretize(ReactionDiffusionSplit(wrapper), FakeDiscretization(), geo)
    @test f.n == 4 * 3
    u₀ = create_initial_condition(f)
    @test length(u₀) == f.n
end

@testset "the trait gate rejects non-participants" begin
    geo = FakeGrid(((0.0, 1.0),), (4,))
    @test_throws ArgumentError semidiscretize(
        ReactionDiffusionSplit(42), FakeDiscretization(), geo,
    )
end
