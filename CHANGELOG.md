# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.0]

Initial extraction of the shared cardiac-EP vocabulary per `DESIGN.md`.

### Added

- `MonodomainModel` — keyword-only, validated continuous declaration with `overrides` (spatial heterogeneity on the model) and `cell_coordinates`.
- `ReactionDiffusionSplit` — pure split annotation with an unconstrained model slot.
- Stimulation vocabulary: `AbstractStimulationProtocol`, abstract `TransmembraneStimulationProtocol`, `NoStimulationProtocol` (strong zero), `AnalyticalTransmembraneStimulationProtocol` (callable + nonzero-interval windows), and the unexported `is_active` query.
- Cell-model contract: `AbstractCellModel`, `is_cell_model`, type-level `num_states` / `state_symbols` / `transmembrane_potential_symbol`, derived `transmembrane_potential_index`, `default_initial_state`, `cell_rhs!`, and the optional `num_parameters` / `parameter_names` / `reaction_rhs!` / `state_rhs!`.
- Capability traits: `has_pointwise_reaction_part`, `reaction_model`, `reaction_solution_symbol`, `reaction_state_symbol`, `reaction_coordinate_system`.
- Pipeline verbs as empty generics with documented contracts: `semidiscretize`, `create_initial_condition`.
