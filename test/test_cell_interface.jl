using CardiacAbstractions
using Test

# ---------------------------------------------------------------------------
# DESIGN.md §3, program 2: a user-defined cell model against the base contract.
# The rhs uses integer literals only, so element-type genericity is real.
# ---------------------------------------------------------------------------
struct MyCell <: AbstractCellModel end
CardiacAbstractions.num_states(::Type{MyCell}) = 2
CardiacAbstractions.state_symbols(::Type{MyCell}) = (:φₘ, :r)
CardiacAbstractions.default_initial_state(::MyCell) = [0.0, 0.0]
function CardiacAbstractions.cell_rhs!(du, u, x, t, m::MyCell)
    # write every slot — never += into an unwritten one
    du[1] = u[1] - u[1]^3 / 3 - u[2]
    du[2] = (u[1] - u[2]) / 5
    return nothing
end

# Only `num_states` implemented: every other query must come from the derived defaults.
struct DefaultCell <: AbstractCellModel end
CardiacAbstractions.num_states(::Type{DefaultCell}) = 3

# Voltage not in slot 1: the derived index must follow the symbols.
struct ReorderedCell <: AbstractCellModel end
CardiacAbstractions.num_states(::Type{ReorderedCell}) = 2
CardiacAbstractions.state_symbols(::Type{ReorderedCell}) = (:r, :φₘ)

# A model that renames the voltage role entirely.
struct RenamedCell <: AbstractCellModel end
CardiacAbstractions.num_states(::Type{RenamedCell}) = 2
CardiacAbstractions.state_symbols(::Type{RenamedCell}) = (:m, :V)
CardiacAbstractions.transmembrane_potential_symbol(::Type{RenamedCell}) = :V

# Voltage symbol absent from the state symbols: the derivation must reject it.
struct BadCell <: AbstractCellModel end
CardiacAbstractions.num_states(::Type{BadCell}) = 2
CardiacAbstractions.state_symbols(::Type{BadCell}) = (:a, :b)

struct Unimplemented <: AbstractCellModel end

struct ParamCell <: AbstractCellModel end
CardiacAbstractions.num_states(::Type{ParamCell}) = 1
CardiacAbstractions.num_parameters(::ParamCell) = 3
CardiacAbstractions.parameter_names(::ParamCell) = (:a, :b, :c)

# Val-wrapping helpers: `Val(f(m))` infers concretely only when `f(m)` folded to a
# compile-time constant — the kernel-facing guarantee the contract makes.
nsval(m) = Val(num_states(m))
idxval(m) = Val(transmembrane_potential_index(m))
symval(m) = Val(state_symbols(m))

@testset "MyCell — §3 program 2 at the declaration level" begin
    m = MyCell()
    @test is_cell_model(m)
    @test num_states(MyCell) == 2
    @test num_states(m) == 2                            # instance → type forwarding
    @test state_symbols(m) == (:φₘ, :r)
    @test transmembrane_potential_symbol(MyCell) == :φₘ  # default fires
    @test transmembrane_potential_index(m) == 1
    @test default_initial_state(m) == [0.0, 0.0]
    model = MonodomainModel(; κ = 1.0e-3, ion = m)       # the program's closing line
    @test model.ion === m
end

@testset "derived defaults" begin
    @test state_symbols(DefaultCell) == (:φₘ, :s1, :s2)
    @test transmembrane_potential_index(DefaultCell()) == 1
    @test transmembrane_potential_index(ReorderedCell()) == 2
    @test transmembrane_potential_symbol(RenamedCell) == :V
    @test transmembrane_potential_index(RenamedCell()) == 2
end

@testset "missing methods error actionably, never StackOverflow" begin
    @test_throws ArgumentError num_states(Unimplemented())
    @test_throws ArgumentError num_states(Unimplemented)   # the type-argument recursion trap
    @test_throws "num_states(::Type{" num_states(Unimplemented)
    # the whole derivation chain bottoms out in the same actionable stub
    @test_throws ArgumentError state_symbols(Unimplemented)
    @test_throws ArgumentError transmembrane_potential_index(Unimplemented())
    du = fill(NaN, 2)
    @test_throws ArgumentError cell_rhs!(du, [0.0, 0.0], nothing, 0.0, Unimplemented())
end

@testset "voltage symbol must appear in state_symbols" begin
    err = try
        transmembrane_potential_index(BadCell())
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin(":φₘ", err.msg)
    @test occursin("state_symbols", err.msg)
end

@testset "type-level queries fold to compile-time constants (no @generated)" begin
    @test @inferred(nsval(MyCell())) === Val(2)
    @test @inferred(idxval(MyCell())) === Val(1)
    @test @inferred(symval(MyCell())) === Val((:φₘ, :r))
    @test @inferred(idxval(DefaultCell())) === Val(1)     # folds through the derived default
    @test @inferred(idxval(ReorderedCell())) === Val(2)
end

@testset "cell_rhs! writes every slot" begin
    m = MyCell()
    u = [0.1, 0.2]
    du = fill(NaN, 2)   # an unwritten slot survives as NaN
    @test cell_rhs!(du, u, nothing, 0.0, m) === nothing
    @test all(!isnan, du)
    @test du[1] ≈ 0.1 - 0.1^3 / 3 - 0.2
    @test du[2] ≈ (0.1 - 0.2) / 5

    # `x` is opaque to a coordinate-free model — any position must be accepted
    du2 = fill(NaN, 2)
    cell_rhs!(du2, u, (0.5, 1.0), 0.0, m)
    @test du2 == du
end

@testset "cell_rhs! computes in eltype(u)" begin
    m = MyCell()
    u64 = [0.1, 0.2]
    du64 = fill(NaN, 2)
    cell_rhs!(du64, u64, nothing, 0.0, m)

    u32 = Float32[0.1, 0.2]
    du32 = fill(NaN32, 2)
    @test cell_rhs!(du32, u32, nothing, 0.0f0, m) === nothing
    @test all(!isnan, du32)
    @test du32 ≈ Float32.(du64)
    @inferred cell_rhs!(du32, u32, nothing, 0.0f0, m)
end

@testset "DiffEq functor convenience — coordinate-free path only" begin
    m = MyCell()
    u = [0.3, -0.1]
    du_direct = fill(NaN, 2)
    cell_rhs!(du_direct, u, nothing, 0.0, m)
    du_functor = fill(NaN, 2)
    @test m(du_functor, u, nothing, 0.0) === nothing
    @test du_functor == du_direct
    # a non-`nothing` p keeps its semantics with the model's owner
    @test_throws MethodError m(du_functor, u, (a = 1,), 0.0)
end

@testset "optional methods stay hasmethod-discoverable" begin
    @test !hasmethod(num_parameters, Tuple{MyCell})
    @test !hasmethod(parameter_names, Tuple{MyCell})
    @test hasmethod(num_parameters, Tuple{ParamCell})
    @test num_parameters(ParamCell()) == 3
    @test parameter_names(ParamCell()) == (:a, :b, :c)
    @test !hasmethod(
        reaction_rhs!,
        Tuple{Vector{Float64}, Float64, Vector{Float64}, Nothing, Float64, MyCell},
    )
    @test !hasmethod(
        state_rhs!,
        Tuple{Vector{Float64}, Float64, Vector{Float64}, Nothing, Float64, MyCell},
    )
end

@testset "is_cell_model" begin
    @test is_cell_model(MyCell())
    @test !is_cell_model(42)
    @test !is_cell_model("not a cell")
end
