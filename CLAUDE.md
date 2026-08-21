# CardiacAbstractions.jl

Zero-dependency base package owning the continuum-level cardiac-EP vocabulary shared by Lightning.jl and Thunderbolt.jl. `DESIGN.md` is the authoritative spec — do not change public shapes without checking it, and do not edit `DESIGN.md` itself without review (it is co-owned with the Thunderbolt maintainer).

## Hard rules

- **Zero dependencies, stdlib only.** A new dependency requires both backend maintainers to agree.
- **Continuum boundary.** Anything touching a mesh, grid, dof, operator, or array layout belongs to a backend, not here.
- **Keyword-only, inner-constructor validation** for multi-slot physics types; no boolean flags in the public surface.
- **Export surface is closed**: the type tree + cell contract + traits + verbs, nothing else (`test/test_exports.jl` pins it). `is_active` is supported API but deliberately unexported.
- Type-level queries must fold to compile-time constants **without `@generated`** (`test_cell_interface.jl` pins folding via the `Val` trick).

## Key files

- `src/cell_interface.jl` — the cell-model contract; the derived voltage index lives here.
- `src/monodomain.jl` — `MonodomainModel` and validation-by-dispatch for the portable field convention.
- `test/test_verbs.jl` — the fake backend proving the §3 canonical program composes.

## Conventions

- Formatting: Runic (CI-enforced, config-free). Tests: SafeTestsets + a JET lint block.
- Test files mirror src files one-to-one; `test_exports.jl` must run before `test_verbs.jl`.
