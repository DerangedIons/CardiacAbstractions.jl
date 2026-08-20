# CardiacAbstractions.jl — Design Specification

**Status:** draft for review · **Date:** 2026-08-20

This document specifies a new base package that becomes the single source of truth for the cardiac-electrophysiology model vocabulary currently duplicated between Lightning.jl (structured grids, matrix-free, GPU) and Thunderbolt.jl (FEM on Ferrite, complex geometry). It specifies the *target* design at the level of abstractions and trade-offs; migration sequencing for either backend, solver design, and kwarg-level surfaces are deliberately out. Thunderbolt's core developer has agreed to depend on such a package, so this doc is written to be reviewable by him per-divergence: every place it departs from a current Thunderbolt shape is flagged inline. Studied inspirations: **Thunderbolt.jl** (donor of most type names, the capability-trait philosophy, and the symbol-derivation design for the voltage role), **Lightning.jl's DESIGN.md** (the congruency-by-convention experiment this package ends, and the keyword-only/validated construction discipline), **CytoZoo.jl** (donor of the cell-model contract's hard-won invariants — the write-every-slot rule, element-type genericity), and **SciMLBase/CommonSolve** (the proven ecosystem pattern for a featherweight package owning shared verbs and declaration types with multiple implementers). Repo home is the DerangedIons org (settled 2026-08-20), co-maintained with the Thunderbolt developer.

## 1. What CardiacAbstractions is

Lightning and Thunderbolt both need `MonodomainModel`, `ReactionDiffusionSplit`, stimulation protocols, a cell-model contract, and a `semidiscretize` verb — and today each owns its own copy, kept congruent by convention. The copies have already drifted: `TransmembraneStimulationProtocol` is abstract in Thunderbolt and concrete in Lightning, the cell-model interfaces (Thunderbolt's `AbstractIonicModel` vs CytoZoo's `AbstractCellModel`) are incompatible, and spatial heterogeneity enters through the split in one package and through the model in the other. CardiacAbstractions is the mesh-free layer both packages import instead: the continuous model declarations, the split annotations, the stimulation vocabulary, the cell-model interface, and the pipeline verbs — with zero dependencies, so neither backend pays anything to take it. The central design claim: the boundary between base and backend is *does it touch a mesh, grid, dof, or array layout* — everything continuous and declarative moves down, everything discrete stays up, and the acceptance test is a user program that switches backends by changing only its geometry and discretization lines.

## 2. Design principles

1. **Declarations, not machinery.** The base owns inert descriptions of the continuous problem; anything that touches a mesh, grid, dof handler, operator, or array layout belongs to a backend.
2. **The contract is methods, not supertypes.** Abstract types exist for convenience and defaults; required paths ask query functions and capability traits, so a type owned by a third package is a full citizen without re-rooting (Thunderbolt's stated philosophy, adopted package-wide).
3. **Zero dependencies.** Julia stdlib only, forever by default; a new dependency requires both backend maintainers to agree, because a dep the base takes is a dep it imposes on everyone.
4. **Portable programs are the acceptance test.** The two-line backend switch in §3 is the contract; a base change that breaks it is wrong even if both backends still compile.
5. **Typed slots, keyword-only, validated at construction.** Physics choices are typed objects in named slots — never flags — constructors are keyword-only (killing Thunderbolt's positional χ/Cₘ swap trap), and validation lives in inner constructors so no path bypasses it.
6. **One name, one meaning.** Every current Lightning/Thunderbolt name clash is resolved here, once; after this doc, a shared word appearing in both backends refers to the base's definition.
7. **Physics-extensible, EP-scoped.** Nothing here designs a second physics (second-implementer rule), but no name or type choice blocks mechanics or perfusion from joining later — which is why the package is named Cardiac, not EP.

## 3. The grammar

```
cell model (CytoZoo | Thunderbolt | yours)  ─┐
stimulation protocol                        ─┼→ MonodomainModel → ReactionDiffusionSplit → semidiscretize(split, disc, geometry) → init / step!
κ, χ, Cₘ  (Number | callable | backend field)┘   [CardiacAbstractions: continuous declaration]      [backend: Lightning | Thunderbolt]    [OrdinaryDiffEqOperatorSplitting]
```

The first canonical program is the reason this package exists — the backend switch. Everything above and below the marked lines is byte-identical:

```julia
using Lightning, CytoZoo                    # or: using Thunderbolt, CytoZoo — both re-export the base vocabulary

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
integrator = init(OperatorSplittingProblem(f, u₀, (0.0, 100.0)),
                  LieTrotterGodunov((Euler(), Euler())); dt = 0.01)
```

The second canonical program is a user-defined cell model — implement the base contract once, run in both backends unchanged:

```julia
struct MyCell <: AbstractCellModel end
CardiacAbstractions.num_states(::Type{MyCell}) = 2
CardiacAbstractions.state_symbols(::Type{MyCell}) = (:φₘ, :r)
CardiacAbstractions.default_initial_state(::MyCell) = [0.0, 0.0]
function CardiacAbstractions.cell_rhs!(du, u, x, t, m::MyCell)
    du[1] = ...; du[2] = ...              # write every slot — never += into an unwritten one
    return nothing
end
model = MonodomainModel(; κ = 1e-3, ion = MyCell())   # works on a grid or a patient mesh
```

The third is spatial heterogeneity in the portable subset — plain values and callables, declared on the model, interpreted by whichever backend discretizes it:

```julia
model = MonodomainModel(; κ = 1e-3, ion = ToRORd(),
    overrides = (celltype = x -> x[1] < 0.5 ? 0.0 : 1.0,))   # endo left, epi right — both backends honor it
```

## 4. The type system

```
AbstractEPModel                                    # convenience family root; the real contract is trait-based
└── MonodomainModel{…}                             #   χ, Cₘ, κ, stim, ion, overrides, cell_coordinates, symbols

ReactionDiffusionSplit{M}                          # pure split annotation; no abstract parent (second-implementer rule)

AbstractStimulationProtocol                        # root of the stimulation vocabulary
└── TransmembraneStimulationProtocol               # abstract: Iₛₜᵢₘ,ᵢ = Iₛₜᵢₘ,ₑ  (resolves the Lightning/Thunderbolt clash)
    ├── NoStimulationProtocol                      #   fieldless — dispatch deletes the term, a strong zero
    └── AnalyticalTransmembraneStimulationProtocol #   plain callable (x,t) + nonzero-interval windows; isbits when f is

AbstractCellModel                                  # convenience root for cell kinetics; contract = the query methods + cell_rhs!

semidiscretize, create_initial_condition           # owned verbs (empty generics); backends add methods
has_pointwise_reaction_part, reaction_* queries    # capability traits semidiscretize asks instead of isa
```

The mockups below are indicative — they fix what each type *owns*, not frozen field inventories; constructors, validation bodies, and defaults live in the implementation.

**`AbstractEPModel` and the capability traits** — the family root exists for dispatch convenience and shared defaults, but the extraction deliberately imports Thunderbolt's deeper rule: code deciding whether something can be split, or what its reaction part is, asks `has_pointwise_reaction_part` / `reaction_model` / `reaction_solution_symbol` / `reaction_state_symbol` / `reaction_coordinate_system` — never `isa AbstractEPModel`. That is what lets a model owned by a third package (an MTK circuit wrapper, a Dict of subdomain models) participate without re-rooting, and it is the seam future physics will reuse: a mechanics model is not an `AbstractEPModel`, but the trait grammar extends without touching this tree.

```julia
abstract type AbstractEPModel end
has_pointwise_reaction_part(model) = false          # capability trait; the EP family answers true
has_pointwise_reaction_part(::AbstractEPModel) = true
reaction_model(model::AbstractEPModel) = model.ion  # + reaction_solution_symbol, reaction_state_symbol,
                                                    #   reaction_coordinate_system — Thunderbolt's names, kept
```

**`MonodomainModel`** — the merged declaration, taking the best half of each parent: Lightning's construction discipline (keyword-only, inner-constructor validation, terse keyword names) and Thunderbolt's field set (including `cell_coordinates`, which a structured-grid backend simply defaults to `nothing`). The one genuine reconciliation: spatial heterogeneity (`overrides`) moves from Lightning's split annotation onto the model, because it is a property of the tissue, not of how the operator is split — Lightning's own DESIGN.md already marks this seam and asks for exactly this move. Coefficient slots stay fully generic; the *portable* value language (§6) is Numbers and callables, and a backend-specific field type in a slot simply pins the model to that backend.

```julia
struct MonodomainModel{Tχ,TC,Tκ,TS<:TransmembraneStimulationProtocol,TI,TO,TX} <: AbstractEPModel
    χ::Tχ                        # surface-to-volume ratio
    Cₘ::TC                       # membrane capacitance per unit area
    κ::Tκ                        # Number = isotropic | NTuple/per-axis | callable x → κ(x) | backend field type
    stim::TS
    ion::TI                      # anything satisfying the cell-model contract below
    overrides::TO                # NamedTuple of spatial parameter fields | nothing (moved here from Lightning's split)
    cell_coordinates::TX         # coordinate the cell model sees as x | nothing = physical (Thunderbolt's slot, kept)
    φ_symbol::Symbol             # published name of the voltage field
    states_symbol::Symbol        # published name of the internal-state block (plural on purpose)
end
```

**`ReactionDiffusionSplit`** — a pure annotation again, now that overrides live on the model: it says *how* to split and nothing else. It exists because the split is a degree of freedom `semidiscretize` cannot infer from the model alone: the same `MonodomainModel` can lower monolithically (one coupled RHS) or as an operator-splitting function whose reaction/diffusion structure the downstream integrator consumes — `semidiscretize(model, disc, geo)` is the monolithic path, and wrapping the model is the declarative request for the split one (Thunderbolt's existing shape, kept). The model parameter is deliberately unconstrained; `semidiscretize` trait-checks `has_pointwise_reaction_part` so foreign model types can be split. No abstract split hierarchy until a second split annotation exists.

```julia
struct ReactionDiffusionSplit{M}
    model::M                     # trait-checked at semidiscretize, not type-constrained here
end
```

**The stimulation vocabulary** — this is where principle 6 bites: Thunderbolt's `TransmembraneStimulationProtocol` is abstract, Lightning's is a concrete callable-plus-windows struct, and the same sentence about "the stimulation protocol" currently means different things in the two codebases. The base keeps Thunderbolt's abstract meaning (the physical statement Iₛₜᵢₘ,ᵢ = Iₛₜᵢₘ,ₑ) and gives Lightning's concrete shape Thunderbolt's concrete name — but with a plain callable where Thunderbolt today has an `AnalyticalCoefficient`, because a coefficient system is discretization machinery: Thunderbolt lowers the callable into its coefficient caches at semidiscretize time, Lightning evaluates it pointwise, and the declaration stays mesh-free. Windows are a tuple of `(t₀, t₁)` pairs rather than a vector of SVectors, which keeps the type isbits and the package dependency-free.

```julia
abstract type AbstractStimulationProtocol end
abstract type TransmembraneStimulationProtocol <: AbstractStimulationProtocol end   # Iₛₜᵢₘ,ᵢ = Iₛₜᵢₘ,ₑ

struct NoStimulationProtocol <: TransmembraneStimulationProtocol end   # fieldless: dispatch deletes the evaluation

struct AnalyticalTransmembraneStimulationProtocol{F,W} <: TransmembraneStimulationProtocol
    f::F                         # callable (x, t) → Iₛₜᵢₘ; positive = depolarizing (the PDE convention, stated once, here)
    nonzero_intervals::W         # NTuple of (t₀, t₁) pairs | nothing; sparsity-in-time hint every backend may exploit
end
```

**The cell-model contract** — the portability keystone, and the piece neither package can supply alone: Thunderbolt's `AbstractIonicModel` and CytoZoo's `AbstractCellModel` are incompatible interfaces to the same concept, which is why a ToRORd built for one backend cannot run in the other today. The base takes CytoZoo's broader root name (`AbstractCellModel` — CytoZoo's zoo includes mitochondria and signaling, not just ion channels) and Thunderbolt's best mechanism: the voltage role is a *symbol* with the index *derived* from `state_symbols`, so names and index cannot disagree, the lookup folds to a compile-time literal, and no `Symbol` survives into a hot kernel. The RHS contract is Thunderbolt's five-argument `cell_rhs!` with the position explicit (and `nothing` allowed for coordinate-free models); heterogeneity does not widen the signature — the backend curries `model.overrides` into the cell model at semidiscretize time (§9 Q1). CytoZoo's DiffEq functor form remains as a one-line convenience wrapper, and CytoZoo's two bitter invariants become part of the documented contract: every `du` slot is written before it is read, and the RHS computes in the element type of `u`. The full required/optional/trait contract routes to a `/julia-interface` follow-up; this doc fixes only the shape:

```julia
abstract type AbstractCellModel end             # convenience root; the contract is the methods, so foreign types
is_cell_model(m) = m isa AbstractCellModel      # opt in via this trait without re-rooting

num_states(::Type{M})                           # required — the first three are properties of the model *type*
state_symbols(::Type{M})                        # default (:φₘ, :s1, :s2, …)
transmembrane_potential_symbol(::Type{M})       # default :φₘ — a role, not a user-facing name
default_initial_state(m)
transmembrane_potential_index(m)                # provided, derived — symbols and index cannot disagree
cell_rhs!(du, u, x, t, m)                       # required; x may be nothing; write every slot, compute in eltype(u)

num_parameters(m), parameter_names(m)           # optional reflection, atomic models only (fitting, parameter studies)
reaction_rhs!(dφ, φ, s, x, t, m)                # optional sub-split forms for backends splitting within the cell
state_rhs!(ds, φ, s, x, t, m)
```

**The verbs** — `semidiscretize(model_or_split, discretization, geometry)` and `create_initial_condition(f)` are owned by the base as empty generic functions with documented contracts, CommonSolve-style: the base fixes the argument grammar and the promise (a spatially-discrete, time-continuous function; an initial state allocated where `f` lives), each backend adds methods for its own discretization and geometry types. This turns today's congruency-by-convention into congruency-by-import — two backends defining methods of the *same function* cannot drift on its meaning.

```julia
function semidiscretize end          # (model | split, discretization, geometry) → semidiscrete function
function create_initial_condition end # (f) → u₀, allocated on the device/layout f lives on
```

## 5. Alternatives considered

**Status quo: congruency by convention.** Zero coupling, total independence, and it survived Thunderbolt's 0.0.x churn precisely because nothing imported anything. It loses because the drift has already happened (the protocol name clash, two cell interfaces), every future Thunderbolt model change needs a hand-mirrored Lightning change, and — decisively — the reason for choosing it evaporated when Thunderbolt's core developer agreed to share a base.

**Lightning depends on Thunderbolt directly.** Single source with no new package. Loses immediately on weight: Thunderbolt carries Ferrite, LinearSolve, JLD2, WriteVTK, and a custom-registry requirement, all of which Lightning's featherweight structured-grid stack exists to avoid; and it makes Lightning hostage to Thunderbolt's release cadence rather than to a tiny stable contract.

**CytoZoo as the shared cell layer, base owns tissue types only.** Thunderbolt takes a CytoZoo dependency and drops its ionic models. Cleanest single source for cells, but it couples the base effort to CytoZoo's coupling-graph machinery and release cadence, and it is a much larger ask of the Thunderbolt developer than depending on a zero-dep interface package. Settled 2026-08-20: the base owns the *interface*; CytoZoo and Thunderbolt both conform, and their concrete models become interchangeable without either depending on the other.

**A traits-and-verbs-only package (no concrete structs).** The smallest possible ask upstream. Loses because the concrete declarations *are* the duplication — `MonodomainModel` is byte-similar in both packages today, it is entirely mesh-free, and leaving it duplicated preserves exactly the drift this package exists to end.

**A symbolic front-end (MTK/SciCompDSL) as the shared language.** Attractive for model exchange in the abstract, but CytoZoo already dropped MTK for cause at the cell level, Lightning's value is hand-tuned matrix-free kernels a symbolic layer cannot see, and Thunderbolt keeps MTK strictly behind extensions. A numerics-first struct vocabulary is the shape both backends actually consume.

## 6. Cross-cutting rules

**Portable field convention.** In any coefficient or override slot: a `Number` means homogeneous, a callable `x -> v` or `(x, t) -> v` means an analytic field, and both are honored by every backend — that is the portability contract. Any other type is backend-specific and pins the model to that backend; backends must reject unknown field types at `semidiscretize`, not at solve time.

**Construction.** Keyword-only public constructors; validation (positivity of χ, Cₘ, κ; symbol distinctness) in inner constructors so generated positional constructors cannot bypass it; no boolean flags anywhere in the public surface.

**Kernel-facing invariants.** `cell_rhs!` writes every slot of `du` before any slot is read-modified, computes in `eltype(u)` (no `Float64` literals leaking into a `Float32` run), and all type-level queries (`num_states`, `state_symbols`, the voltage role) fold to compile-time constants — deliberately without `@generated`, for the world-age reason documented in Thunderbolt.

**Sign and units.** Stimulus positive = depolarizing, stated once here and inherited by both backends; no unit system is fixed — κ, χ, Cₘ, the stimulus, and the cell model must simply agree, and the base documents (as Lightning does) that κ enters as κ/(χCₘ) and the stimulus as Iₛₜᵢₘ/Cₘ.

**Versioning and distribution.** Registered in **General** (Thunderbolt lives in a custom registry, but a General-registered base is reachable by both and keeps Lightning General-only). The base is the compatibility contract between backends: breaking releases are rare, deliberate, and coordinated so both backends cross a major together. Both backends `Reexport` the base vocabulary — end users write `using Lightning` or `using Thunderbolt` and never import CardiacAbstractions directly.

**Export surface.** Everything in §4's tree plus the verbs and traits; nothing else. The base never exports a solver, a mesh, a layout, or a concrete cell model.

## 7. Non-goals

**A second physics (mechanics, perfusion, fluids).** The trait grammar and the Cardiac-not-EP name keep the door open; designing mechanics abstractions before a second implementer exists violates the second-implementer rule. Thunderbolt's internal-variable machinery stays put until then.

**Bidomain declarations** (parabolic–elliptic and parabolic–parabolic). Real vocabulary — Thunderbolt carries `::Any`-field placeholder structs today — but no backend can discretize either form, so the properly-typed declarations join the base when one can (the same second-implementer rule). `AbstractStimulationProtocol` keeps its wider root as the seam split intra-/extracellular stimulation re-enters through. Dropped from v1 in review, 2026-08-20.

**Concrete cell models.** None in the base, ever — CytoZoo and Thunderbolt keep their zoos; the shared interface makes them interchangeable, which lets today's duplicated FHN/AlievPanfilov wither naturally (§9 Q3).

**Solvers, problems, integrators.** OrdinaryDiffEqOperatorSplitting is already the shared downstream driver for both backends; the base owns nothing OS owns.

**A coefficient system.** The portable convention (§6) is the whole story; Thunderbolt's coefficient caches and Lightning's spatial functors are backend machinery.

**Solution access vocabulary** (`getvariable`, SymbolicIndexingInterface glue). Phase 2: the model symbols standardized here are the prerequisite; the accessor verbs migrate once both backends' solution shapes stabilize.

**ECG and postprocessing.** Lead-field computation needs a mesh; stays in Thunderbolt. A declaration-level ECG model may migrate later.

**Rush–Larsen gate structure.** Both backends want RL, but exposing gate/α/β structure through the base interface is deferred until both RL paths are stable enough to name a shared contract.

**Coordinate systems and microstructure** (heart axes, LV coordinates, fiber fields). Geometric, therefore backend; `cell_coordinates` on the model is the declaration-level hook and is all the base knows.

## 8. Ecosystem alignment

The package follows the SciMLBase/CommonSolve school: a near-zero-dep base owning declaration types and verb contracts, with heavyweight implementers adding methods — the one pattern in the Julia ecosystem with a decade of evidence that independent packages can share a vocabulary without sharing a release cadence. Within the cardiac stack it deliberately inverts Thunderbolt's current layering in one respect: Thunderbolt's stimulation protocols root under its `AbstractSourceTerm` and wrap coefficient types, which the base cannot know about — so Thunderbolt will need to treat source-term-ness as a trait or interface rather than a supertype for the shared protocol family (§9 Q2, the one migration item that needs the core developer's explicit sign-off). CytoZoo re-roots its `AbstractCellModel` under the base and gains Thunderbolt compatibility for its entire zoo in return; Lightning's planned `StateBlockedLayout`, recorders, and MFO glue are untouched — they sit exactly on the backend side of the boundary this doc draws.

## 9. Open questions

1. **Heterogeneity seam in `cell_rhs!`.** (a) Keep the five-argument signature and have backends curry `model.overrides` into the cell model at semidiscretize time via an optional hook (e.g. `bind_overrides(ion, overrides)`, error by default); (b) widen the signature with an explicit parameter-context argument `cell_rhs!(du, u, p, x, t, m)`. **Recommend (a)** — the hot signature stays minimal, homogeneous models pay nothing, and CytoZoo's SpatialContext becomes an implementation detail of its bridge; the full contract gets its `/julia-interface` pass before implementation.
2. **Thunderbolt's `AbstractSourceTerm` rooting.** The base's protocol family cannot subtype a Thunderbolt-internal abstract type. (a) Thunderbolt converts source-term-ness to a trait; (b) `AbstractSourceTerm` itself moves into the base. **Recommend (a)** — moving it drags element-cache semantics into a package that must stay mesh-free — but this is the core developer's call.
3. **Fate of the duplicated toy cell models** (FHN, AlievPanfilov in both Thunderbolt and CytoZoo). (a) Leave both, interchangeable under the shared interface, and let usage decide; (b) actively consolidate into CytoZoo. **Recommend (a)** — consolidation is a governance negotiation with no technical payoff once the interface makes them equivalent.
4. **Verb name `create_initial_condition`.** (a) Adopt Lightning's name as the base verb; (b) something shorter (`initial_state`) at the cost of renaming in both backends. **Recommend (a)** — it is already in the wild in Lightning and unclaimed in Thunderbolt.
