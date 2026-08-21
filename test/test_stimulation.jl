using CardiacAbstractions
using Test

const pulse = (x, t) -> (t <= 2.0 && x[1] <= 1.5) ? 50.0 : 0.0

function windowed_protocol()
    return AnalyticalTransmembraneStimulationProtocol(;
        f = pulse,
        nonzero_intervals = ((0.0, 2.0), (500.0, 502.0)),
    )
end

unwindowed_protocol() = AnalyticalTransmembraneStimulationProtocol(; f = pulse)

struct ForeignProtocol <: AbstractStimulationProtocol end

@testset "hierarchy resolves the protocol name clash" begin
    @test isabstracttype(TransmembraneStimulationProtocol)
    @test TransmembraneStimulationProtocol <: AbstractStimulationProtocol
    @test AnalyticalTransmembraneStimulationProtocol <: TransmembraneStimulationProtocol
    @test NoStimulationProtocol <: TransmembraneStimulationProtocol
end

@testset "keyword-only construction" begin
    p = unwindowed_protocol()
    @test p.f === pulse
    @test p.nonzero_intervals === nothing
    @test_throws MethodError AnalyticalTransmembraneStimulationProtocol(pulse)
    @test_throws UndefKeywordError AnalyticalTransmembraneStimulationProtocol()
end

@testset "interval normalization and rejection" begin
    @test windowed_protocol().nonzero_intervals === ((0.0, 2.0), (500.0, 502.0))

    mixed = AnalyticalTransmembraneStimulationProtocol(;
        f = pulse,
        nonzero_intervals = ((0, 2.0), (5, 6)),
    )
    @test mixed.nonzero_intervals == ((0, 2.0), (5, 6))

    from_vector = AnalyticalTransmembraneStimulationProtocol(;
        f = pulse,
        nonzero_intervals = [(0.0, 2.0)],
    )
    @test from_vector.nonzero_intervals isa Tuple
    @test from_vector.nonzero_intervals == ((0.0, 2.0),)

    degenerate = AnalyticalTransmembraneStimulationProtocol(;
        f = pulse,
        nonzero_intervals = ((1.0, 1.0),),
    )
    @test CardiacAbstractions.is_active(degenerate, 1.0)

    @test_throws ArgumentError AnalyticalTransmembraneStimulationProtocol(;
        f = pulse,
        nonzero_intervals = ((2.0, 0.0),),
    )
    @test_throws ArgumentError AnalyticalTransmembraneStimulationProtocol(;
        f = pulse,
        nonzero_intervals = (),
    )
    @test_throws ArgumentError AnalyticalTransmembraneStimulationProtocol(;
        f = pulse,
        nonzero_intervals = ((0.0, 1.0, 2.0),),
    )
end

@testset "isbits exactly when f is" begin
    @test isbitstype(typeof(windowed_protocol()))
    @test isbitstype(NoStimulationProtocol)
    @test Base.issingletontype(NoStimulationProtocol)
    captured = rand(3)
    closure_protocol = AnalyticalTransmembraneStimulationProtocol(; f = (x, t) -> captured[1])
    @test !isbitstype(typeof(closure_protocol))
end

@testset "functor evaluation and the strong zero" begin
    p = windowed_protocol()
    x = (1.0, 0.0, 0.0)
    @test p(x, 1.0) == 50.0
    @test p(x, 3.0) == 0.0
    @test p((2.0, 0.0, 0.0), 1.0) == 0.0
    none = NoStimulationProtocol()
    @test none(x, 1.0) === false           # a strong zero, not an evaluated 0.0
    @test none(x, 1.0) * 123.4 == 0.0
end

@testset "is_active — the sparsity-in-time hint" begin
    is_active = CardiacAbstractions.is_active
    @test !is_active(NoStimulationProtocol(), 0.0)
    always = unwindowed_protocol()
    @test is_active(always, -1.0e6)
    @test is_active(always, 1.0e6)
    w = windowed_protocol()
    @test is_active(w, 0.0)                # window endpoints are active
    @test is_active(w, 1.0)
    @test is_active(w, 2.0)
    @test !is_active(w, 2.5)
    @test is_active(w, 501.0)
    @test !is_active(w, 600.0)
    @test is_active(ForeignProtocol(), 0.0)   # conservative fallback for foreign protocols
end

@testset "construction is type-stable" begin
    @test @inferred(windowed_protocol()) isa AnalyticalTransmembraneStimulationProtocol
    @test @inferred(unwindowed_protocol()) isa AnalyticalTransmembraneStimulationProtocol
end
