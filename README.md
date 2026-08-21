<h1 align="center">CardiacAbstractions.jl</h1>

<p align="center"><em>The mesh-free vocabulary of cardiac electrophysiology — one set of model declarations, many backends.</em></p>

<p align="center">
  <a href="https://DerangedIons.github.io/CardiacAbstractions.jl/stable"><img src="https://img.shields.io/badge/docs-stable-blue.svg" alt="Stable Docs"></a>
  <a href="https://DerangedIons.github.io/CardiacAbstractions.jl/dev"><img src="https://img.shields.io/badge/docs-dev-blue.svg" alt="Dev Docs"></a>
  <a href="https://github.com/DerangedIons/CardiacAbstractions.jl/actions/workflows/CI.yml?query=branch%3Amain"><img src="https://github.com/DerangedIons/CardiacAbstractions.jl/actions/workflows/CI.yml/badge.svg?branch=main" alt="Build Status"></a>
  <a href="https://codecov.io/gh/DerangedIons/CardiacAbstractions.jl"><img src="https://codecov.io/gh/DerangedIons/CardiacAbstractions.jl/graph/badge.svg" alt="codecov"></a>
  <a href="https://github.com/DerangedIons/CardiacAbstractions.jl/blob/main/LICENSE"><img src="https://img.shields.io/badge/license-MIT-blue" alt="License"></a>
</p>

---

## Why CardiacAbstractions?

Cardiac-electrophysiology backends keep reinventing the same vocabulary: [Lightning.jl](https://github.com/DerangedIons/Lightning.jl) (structured grids, matrix-free, GPU) and [Thunderbolt.jl](https://github.com/termi-official/Thunderbolt.jl) (FEM on Ferrite, complex geometry) each carried their own `MonodomainModel`, stimulation protocols, cell-model interface, and `semidiscretize` verb — kept congruent only by convention, and already drifting. CardiacAbstractions.jl is the zero-dependency base both import instead: the continuous model declarations, the split annotations, the stimulation vocabulary, the cell-model contract, and the pipeline verbs as empty generic functions each backend adds methods to, CommonSolve-style.

The boundary rule: anything touching a mesh, grid, dof handler, operator, or array layout belongs to a backend; everything continuous and declarative lives here. The acceptance test is a user program that switches backends by changing only its geometry and discretization lines.

## Quick Start

```julia
using Lightning, CytoZoo                    # or: using Thunderbolt, CytoZoo — both re-export this vocabulary

model = MonodomainModel(; κ = (0.133, 0.0176, 0.0176), χ = 140.0, Cₘ = 0.01,
    ion  = TenTusscher2006(),
    stim = AnalyticalTransmembraneStimulationProtocol(; f = (x, t) -> (t ≤ 2.0 && x[1] ≤ 1.5) ? 50.0 : 0.0,
                                                        nonzero_intervals = ((0.0, 2.0),)))
split = ReactionDiffusionSplit(model)

geo  = CartesianGrid(((0.0, 20.0), (0.0, 7.0), (0.0, 3.0)), (100, 35, 15))   # ← backend line 1
f    = semidiscretize(split, FiniteDifferenceDiscretization(), geo)          # ← backend line 2
# Thunderbolt instead:  geo = load_mesh("patient_lv.msh")
#                       f   = semidiscretize(split, FiniteElementDiscretization(...), geo)

u₀ = create_initial_condition(f)
```

Or implement the cell-model contract once and run your model in both backends unchanged:

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

## What's in the box

- `MonodomainModel` — the merged continuous declaration: keyword-only, validated at construction, with spatial heterogeneity (`overrides`) declared on the model.
- `ReactionDiffusionSplit` — a pure split annotation; `semidiscretize` trait-checks the model, so foreign model types can be split.
- The stimulation vocabulary — abstract `TransmembraneStimulationProtocol`, the fieldless `NoStimulationProtocol` strong zero, and `AnalyticalTransmembraneStimulationProtocol` (a plain callable plus optional nonzero-interval windows; isbits when the callable is).
- The cell-model contract — `AbstractCellModel`, type-level queries that fold to compile-time constants, and a symbol-derived voltage index so names and index cannot disagree.
- The verbs — `semidiscretize` and `create_initial_condition` as empty generics with documented contracts.

What's deliberately **not** here: solvers, meshes, layouts, coefficient systems, and concrete cell models — those belong to the backends and to CytoZoo.

## Installation

```julia
using Pkg
Pkg.add("CardiacAbstractions")
```

Requires Julia 1.10 or newer. End users typically never install this directly — Lightning and Thunderbolt re-export the vocabulary.

## Documentation

- [Getting started and the grammar](https://DerangedIons.github.io/CardiacAbstractions.jl/dev)
- [The cell-model contract](https://DerangedIons.github.io/CardiacAbstractions.jl/dev/cell_models)
- [API reference](https://DerangedIons.github.io/CardiacAbstractions.jl/dev/api)
- `DESIGN.md` — the full design specification: boundary, trade-offs, alternatives considered.

## Contributing

Issues and pull requests are welcome. This package is the compatibility contract between backends, so breaking changes are rare, deliberate, and coordinated with both backend maintainers.
