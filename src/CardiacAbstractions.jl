"""
    CardiacAbstractions

The continuum-level vocabulary shared by cardiac-electrophysiology backends: continuous model
declarations, split annotations, the stimulation vocabulary, the cell-model contract, and
the pipeline verbs. The boundary rule is simple — anything touching a mesh, grid, dof
handler, operator, or array layout belongs to a backend (Lightning, Thunderbolt); the
inert description of the continuous problem lives here, so a user program switches
backends by changing only its geometry and discretization lines.

Zero dependencies, forever by default: a dependency the base takes is a dependency it
imposes on every backend.
"""
module CardiacAbstractions

include("cell_interface.jl")
export AbstractCellModel, is_cell_model
export num_states, state_symbols, transmembrane_potential_symbol,
    transmembrane_potential_index, default_initial_state, cell_rhs!
export num_parameters, parameter_names, reaction_rhs!, state_rhs!

include("stimulation.jl")
export AbstractStimulationProtocol, TransmembraneStimulationProtocol,
    NoStimulationProtocol, AnalyticalTransmembraneStimulationProtocol

include("traits.jl")
export AbstractEPModel
export has_pointwise_reaction_part, reaction_model, reaction_solution_symbol,
    reaction_state_symbol, reaction_coordinate_system

include("monodomain.jl")
export MonodomainModel

include("split.jl")
export ReactionDiffusionSplit

include("verbs.jl")
export semidiscretize, create_initial_condition

end # module
