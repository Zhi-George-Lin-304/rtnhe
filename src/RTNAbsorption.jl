module RTNAbsorption

using LinearAlgebra
using DifferentialEquations
using Plots
using Printf
using Random

include("utils.jl")
include("spectrum.jl")
include("rtn_hierarchy.jl")
include("gaussian_benchmark.jl")
include("channels.jl")
include("models.jl")
include("simulation.jl")
include("plotting.jl")
include("examples.jl")

export
    logspace_int,
    trapezoid_weights,
    safe_normalize!,
    effective_channel_number,
    total_noise_variance,

    absorption_spectrum,

    RTNGroupedCoherenceParams,
    make_rtn_grouped_coherence_params,
    rtn_grouped_coherence_rhs!,
    solve_rtn_absorption_response,

    gaussian_factor_grouped,
    gaussian_absorption_response,

    _sample_categorical,
    normalize_rtn_deltas,
    sample_rtn_counts,
    sample_rtn_gammas,
    sample_rtn_amplitude_weights,

    identical_channels,
    typed_distinct_groups,
    unique_distinct_channels,

    default_single_transition_model,
    solve_and_spectrum,

    make_identical_absorption_plot,
    make_identical_absorption_error_plot,
    make_identical_absorption_gif,
    make_all_identical_absorption_outputs,

    make_typed_distinct_absorption_plot,
    make_typed_distinct_absorption_error_plot,
    make_typed_distinct_absorption_gif,
    make_all_typed_distinct_absorption_outputs,

    make_unique_distinct_absorption_plot,
    make_unique_distinct_absorption_error_plot,
    make_unique_distinct_absorption_gif,
    make_all_unique_distinct_absorption_outputs,

    run_identical_examples,
    run_typed_examples,
    run_typed_random_examples,
    run_typed_weighted_examples,
    run_unique_examples,
    run_unique_random_examples,
    main

end