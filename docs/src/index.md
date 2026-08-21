# CardiacAbstractions.jl

CardiacAbstractions.jl is the continuum-level vocabulary shared by cardiac-electrophysiology backends: the continuous model declarations (`MonodomainModel`), the split annotations (`ReactionDiffusionSplit`), the stimulation vocabulary, the cell-model contract (`AbstractCellModel` and its query methods), and the pipeline verbs (`semidiscretize`, `create_initial_condition`) as empty generic functions each backend adds methods to.

The boundary rule is simple: **anything touching a mesh, grid, dof handler, operator, or array layout belongs to a backend** — Lightning.jl (structured grids, matrix-free, GPU) or Thunderbolt.jl (FEM on Ferrite, complex geometry). Everything continuous and declarative lives here, with zero dependencies, so neither backend pays anything to take it. Both backends re-export this vocabulary, so end users write `using Lightning` or `using Thunderbolt` and never import CardiacAbstractions directly.

## The grammar

```
cell model (CytoZoo | Thunderbolt | yours)  ─┐
stimulation protocol                        ─┼→ MonodomainModel → ReactionDiffusionSplit → semidiscretize(split, disc, geometry) → init / step!
κ, χ, Cₘ  (Number | callable | backend field)┘   [CardiacAbstractions: continuous declaration]      [backend: Lightning | Thunderbolt]    [OrdinaryDiffEqOperatorSplitting]
```

## The backend switch

The acceptance test of the whole design: a user program switches backends by changing only its geometry and discretization lines. Everything above and below the marked lines is byte-identical.

```julia
using Lightning, CytoZoo                    # or: using Thunderbolt, CytoZoo

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

## Spatial heterogeneity

Heterogeneity is declared on the model — it is a property of the tissue, not of how the operator is split. The portable value language is `Number`s and callables; both are honored by every backend, and any other type pins the model to the backend that understands it.

```julia
model = MonodomainModel(; κ = 1.0e-3, ion = ToRORd(),
    overrides = (celltype = x -> x[1] < 0.5 ? 0.0 : 1.0,))   # endo left, epi right
```

## Sign and units

The stimulus is positive = depolarizing, stated once in this package and inherited by every backend. No unit system is fixed — κ, χ, Cₘ, the stimulus, and the cell model must simply agree; κ enters the equation as κ/(χCₘ) and the stimulus as Iₛₜᵢₘ/Cₘ.

## What lives elsewhere

No solvers, meshes, layouts, coefficient systems, or concrete cell models — ever. See the design specification (`DESIGN.md` in the repository) for the full boundary, the alternatives considered, and the open questions.
