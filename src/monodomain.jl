# Validation by dispatch implements the portable field convention: values the base can
# check are checked; values only a backend can interpret (callables, backend field types)
# pass through, and backends reject unknown field types at `semidiscretize`, not at solve
# time.
function _validate_coefficient(name, v::Number)
    v > zero(v) || throw(ArgumentError("$name must be positive, got $v"))
    return nothing
end

function _validate_coefficient(name, v::Union{Tuple{Vararg{Number}}, AbstractVector{<:Number}})
    isempty(v) && throw(ArgumentError("per-axis $name must not be empty"))
    all(x -> x > zero(x), v) || throw(
        ArgumentError("every per-axis $name entry must be positive, got $v")
    )
    return nothing
end

_validate_coefficient(name, v) = nothing

"""
    MonodomainModel(; κ, ion, χ = 1, Cₘ = 1, stim = NoStimulationProtocol(),
                      overrides = nothing, cell_coordinates = nothing,
                      φ_symbol = :φₘ, states_symbol = :states)

The monodomain equation

```
χCₘ∂ₜφₘ = ∇⋅κ∇φₘ + χ(Iᵢₒₙ(φₘ, s, t) + Iₛₜᵢₘ(x, t))
                 ∂ₜs = f(φₘ, s, t)
```

for transmembrane potential `φₘ` and cell-model states `s`, declared at the continuum
level. A model
is an inert description of the continuous problem; [`semidiscretize`](@ref) turns it into
a right-hand side on whatever geometry a backend provides.

# Keyword Arguments

  - `κ`: conductivity. A `Number` is isotropic; an `NTuple`/`AbstractVector` is per-axis
    anisotropy; a callable `x -> κ(x)` is an analytic field. Numbers and callables are the
    portable subset every backend honors — any other type is backend-specific and pins the
    model to that backend, which rejects unknown field types at `semidiscretize`.
  - `ion`: anything satisfying the cell-model contract (see [`AbstractCellModel`](@ref)).
    Supplies `Iᵢₒₙ` and the state kinetics; it *is* the reaction part.
  - `χ`: surface-to-volume ratio.
  - `Cₘ`: membrane capacitance per unit area.
  - `stim`: a [`TransmembraneStimulationProtocol`](@ref). Positive = depolarizing.
  - `overrides`: a `NamedTuple` of spatial parameter fields the cell model resolves
    pointwise (cell type, pH, hypoxia, …), or `nothing`. Values follow the portable
    convention: `Number`, callable, or a backend field type. Heterogeneity lives on the
    model — it is a property of the tissue, not of how the operator is split — and
    backends curry it into the cell model at `semidiscretize` time.
  - `cell_coordinates`: the coordinate the cell model sees as `x`, or `nothing` for
    physical coordinates. Interpreting a non-`nothing` value is the backend's job.
  - `φ_symbol`, `states_symbol`: published names of the voltage field and of the
    internal-state block, surfaced through [`reaction_solution_symbol`](@ref) and
    [`reaction_state_symbol`](@ref). `states_symbol` is deliberately plural: a cell model
    may itself name a state `:s`, and a singular default would shadow it.

The constructor is **keyword-only on purpose**: `χ` and `Cₘ` are dimensionally
interchangeable in the equation above, so a swapped positional pair would be silent.
There is no positional constructor at all — validation lives in the inner constructor and
no path bypasses it.

# Units

Nothing here fixes a unit system — `κ`, `χ`, `Cₘ`, the stimulus, and the cell model must
simply agree. What is fixed is that `κ` enters the equation as the diffusivity `κ/(χCₘ)`
and the stimulus as `Iₛₜᵢₘ/Cₘ`.

### Examples

```julia
model = MonodomainModel(; κ = 1.0e-3, ion = MyCell())

# fibres along x, three-dimensional anisotropy
model = MonodomainModel(;
    κ = (0.133, 0.0176, 0.0176),
    χ = 140.0,
    Cₘ = 0.01,
    ion = MyCell(),
)

# spatial heterogeneity in the portable subset: endo left, epi right
model = MonodomainModel(;
    κ = 1.0e-3,
    ion = MyCell(),
    overrides = (celltype = x -> x[1] < 0.5 ? 0.0 : 1.0,),
)
```

See also: [`ReactionDiffusionSplit`](@ref), [`semidiscretize`](@ref).
"""
struct MonodomainModel{Tχ, TC, Tκ, TS <: TransmembraneStimulationProtocol, TI, TO, TX} <: AbstractEPModel
    χ::Tχ
    Cₘ::TC
    κ::Tκ
    stim::TS
    ion::TI
    overrides::TO
    cell_coordinates::TX
    φ_symbol::Symbol
    states_symbol::Symbol

    # Keyword-only inner constructor: defining it suppresses the generated positional
    # constructors, so no construction path bypasses validation.
    function MonodomainModel(;
            κ,
            ion,
            χ = 1,
            Cₘ = 1,
            stim::TransmembraneStimulationProtocol = NoStimulationProtocol(),
            overrides::Union{NamedTuple, Nothing} = nothing,
            cell_coordinates = nothing,
            φ_symbol::Symbol = :φₘ,
            states_symbol::Symbol = :states,
        )
        _validate_coefficient("conductivity κ", κ)
        _validate_coefficient("χ", χ)
        _validate_coefficient("Cₘ", Cₘ)
        φ_symbol === states_symbol && throw(
            ArgumentError("φ_symbol and states_symbol must differ, both are :$φ_symbol")
        )
        return new{
            typeof(χ), typeof(Cₘ), typeof(κ), typeof(stim), typeof(ion),
            typeof(overrides), typeof(cell_coordinates),
        }(
            χ, Cₘ, κ, stim, ion, overrides, cell_coordinates, φ_symbol, states_symbol
        )
    end
end
