using CardiacAbstractions
using Test

# No concrete cell models in the base, ever — a stub stands in for CytoZoo's zoo.
struct TenTusscher2006 <: AbstractCellModel end
CardiacAbstractions.num_states(::Type{TenTusscher2006}) = 2
CardiacAbstractions.default_initial_state(::TenTusscher2006) = [0.0, 1.0]
function CardiacAbstractions.cell_rhs!(du, u, x, t, m::TenTusscher2006)
    du[1] = -u[1]
    du[2] = -u[2]
    return nothing
end

struct ForeignProtocol <: AbstractStimulationProtocol end   # not transmembrane
struct FakeBackendField end

function canonical_model()
    return MonodomainModel(;
        κ = (0.133, 0.0176, 0.0176),
        χ = 140.0,
        Cₘ = 0.01,
        ion = TenTusscher2006(),
        stim = AnalyticalTransmembraneStimulationProtocol(;
            f = (x, t) -> (t ≤ 2.0 && x[1] ≤ 1.5) ? 50.0 : 0.0,
            nonzero_intervals = ((0.0, 2.0),),
        ),
    )
end

@testset "§3 program 1 — the declaration half" begin
    model = canonical_model()
    @test model isa MonodomainModel
    @test model.κ === (0.133, 0.0176, 0.0176)
    @test model.χ === 140.0
    @test model.Cₘ === 0.01
    @test model.ion isa TenTusscher2006
    @test model.stim isa AnalyticalTransmembraneStimulationProtocol
    @test model.stim((1.0, 0.0, 0.0), 1.0) == 50.0
    @test model.stim((1.0, 0.0, 0.0), 3.0) == 0.0
    split = ReactionDiffusionSplit(model)
    @test split.model === model
end

@testset "keyword-only: no positional path exists" begin
    ion = TenTusscher2006()
    @test_throws MethodError MonodomainModel(
        140.0, 0.01, 1.0e-3, NoStimulationProtocol(), ion, nothing, nothing, :φₘ, :states,
    )
    @test_throws UndefKeywordError MonodomainModel(; κ = 1.0e-3)    # ion is required
    @test_throws UndefKeywordError MonodomainModel(; ion = ion)     # κ is required
end

@testset "rejected configurations" begin
    ion = TenTusscher2006()
    @test_throws ArgumentError MonodomainModel(; κ = 0.0, ion)
    @test_throws ArgumentError MonodomainModel(; κ = -1.0, ion)
    @test_throws ArgumentError MonodomainModel(; κ = (), ion)
    @test_throws ArgumentError MonodomainModel(; κ = (0.1, -0.2), ion)
    @test_throws ArgumentError MonodomainModel(; κ = [0.1, 0.0], ion)
    @test_throws ArgumentError MonodomainModel(; κ = 1.0, χ = 0, ion)
    @test_throws ArgumentError MonodomainModel(; κ = 1.0, χ = -1.0, ion)
    @test_throws ArgumentError MonodomainModel(; κ = 1.0, Cₘ = 0.0, ion)
    @test_throws ArgumentError MonodomainModel(; κ = 1.0, ion, φ_symbol = :u, states_symbol = :u)
    @test_throws TypeError MonodomainModel(; κ = 1.0, ion, overrides = Dict(:a => 1))
    @test_throws TypeError MonodomainModel(; κ = 1.0, ion, stim = "zap")
    # a protocol outside the transmembrane family is rejected by the slot's type
    @test_throws TypeError MonodomainModel(; κ = 1.0, ion, stim = ForeignProtocol())
end

@testset "portable-convention acceptance" begin
    ion = TenTusscher2006()
    @test MonodomainModel(; κ = 1.0e-3, ion) isa MonodomainModel             # Number
    @test MonodomainModel(; κ = (0.1, 0.2), ion) isa MonodomainModel         # per-axis tuple
    @test MonodomainModel(; κ = [0.1, 0.2], ion) isa MonodomainModel         # per-axis vector
    @test MonodomainModel(; κ = x -> 1.0e-3, ion) isa MonodomainModel        # analytic field
    # backend field types pass through here; the backend rejects unknowns at semidiscretize
    @test MonodomainModel(; κ = FakeBackendField(), ion) isa MonodomainModel
    @test MonodomainModel(; κ = 1.0, χ = x -> 140.0, ion) isa MonodomainModel
    @test MonodomainModel(; κ = 1.0, Cₘ = (x, t) -> 0.01, ion) isa MonodomainModel
end

@testset "§3 program 3 — heterogeneity as overrides on the model" begin
    model = MonodomainModel(;
        κ = 1.0e-3,
        ion = TenTusscher2006(),
        overrides = (celltype = x -> x[1] < 0.5 ? 0.0 : 1.0,),
    )
    @test model.overrides.celltype((0.25,)) == 0.0   # endo left
    @test model.overrides.celltype((0.75,)) == 1.0   # epi right
end

@testset "defaults" begin
    model = MonodomainModel(; κ = 1.0e-3, ion = TenTusscher2006())
    @test model.χ === 1
    @test model.Cₘ === 1
    @test model.stim === NoStimulationProtocol()
    @test model.overrides === nothing
    @test model.cell_coordinates === nothing
    @test model.φ_symbol === :φₘ
    @test model.states_symbol === :states
end

@testset "construction is type-stable" begin
    @test @inferred(canonical_model()) isa MonodomainModel
end
