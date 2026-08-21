# The Cell-Model Contract

The cell-model interface is the portability keystone: implement it once and the model runs in every backend unchanged — on a structured grid or a patient mesh. The contract is the *methods*, not the supertype: `AbstractCellModel` exists for convenience, and a type that already has a supertype in a foreign package opts in via `is_cell_model` instead.

## A complete implementation

```julia
struct MyCell <: AbstractCellModel end
CardiacAbstractions.num_states(::Type{MyCell}) = 2
CardiacAbstractions.state_symbols(::Type{MyCell}) = (:φₘ, :r)
CardiacAbstractions.default_initial_state(::MyCell) = [0.0, 0.0]
function CardiacAbstractions.cell_rhs!(du, u, x, t, m::MyCell)
    du[1] = u[1] - u[1]^3 / 3 - u[2]      # write every slot — never += into an unwritten one
    du[2] = (u[1] - u[2]) / 5
    return nothing
end

model = MonodomainModel(; κ = 1.0e-3, ion = MyCell())   # works on a grid or a patient mesh
```

## Required methods

| Method | Meaning |
|---|---|
| `num_states(::Type{M})` | number of ODE states, ≥ 1; a property of the *type* |
| `default_initial_state(m)` | initial condition, length `num_states(M)` |
| `cell_rhs!(du, u, x, t, m)` | the kinetics; `x` may be `nothing` |

## Optional methods with derived defaults

| Method | Default |
|---|---|
| `state_symbols(::Type{M})` | `(:φₘ, :s1, :s2, …)` |
| `transmembrane_potential_symbol(::Type{M})` | `:φₘ` — a role marker, not a user-facing name |

The voltage role is a *symbol*, and `transmembrane_potential_index` is **derived** from the two queries above — never implemented directly — so names and index cannot disagree. The lookup folds to a compile-time literal, and no `Symbol` survives into a hot kernel.

## Optional methods with no default

Probe these with `hasmethod`; a model that leaves them undefined correctly reports the missing capability.

| Method | Meaning |
|---|---|
| `num_parameters(m)`, `parameter_names(m)` | reflection for atomic models (fitting, parameter studies) |
| `reaction_rhs!(dφ, φ, s, x, t, m)` | sub-split form: the voltage equation alone |
| `state_rhs!(ds, φ, s, x, t, m)` | sub-split form: the non-voltage states alone |

## The two bitter invariants

Both are hard-won and non-negotiable for kernel correctness:

1. **Write every slot of `du` before any slot is read-modified.** `du` arrives uninitialized, so `du[i] += v` as a slot's first write reads garbage — and under a contributory coupling it becomes a silent double-count rather than an obvious `NaN`.
2. **Compute in `eltype(u)`.** Wrap every numeric literal as `T(0.5)` (or use integer literals, which promote losslessly) so no `Float64` intermediate leaks into a `Float32` run. This keeps the same code CPU-, GPU-, and precision-generic.

## Heterogeneity and the stimulus

Spatial parameter heterogeneity never widens the `cell_rhs!` signature: backends curry the tissue model's `overrides` into the cell model at `semidiscretize` time. The tissue model owns the stimulus — a cell model used on the tissue path should carry none of its own (the sign conventions differ, and the two would silently add).

## Foreign types

A type that cannot subtype `AbstractCellModel` opts in without re-rooting:

```julia
struct ForeignIon end   # supertype owned by another package
CardiacAbstractions.is_cell_model(::ForeignIon) = true
CardiacAbstractions.num_states(::Type{ForeignIon}) = 2
CardiacAbstractions.default_initial_state(::ForeignIon) = [0.0, 0.0]
CardiacAbstractions.cell_rhs!(du, u, x, t, m::ForeignIon) = ...
```
