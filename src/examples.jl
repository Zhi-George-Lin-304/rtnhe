include("plotting.jl")

# Each mode below produces three outputs:
#   1. PNG spectrum plot
#   2. PNG error-vs-n plot
#   3. GIF spectrum animation

const DEFAULT_MAX_HIERARCHY_SIZE = 2_000_000


function run_identical_examples(; max_hierarchy_size = DEFAULT_MAX_HIERARCHY_SIZE)
    make_identical_absorption_plot(
        ns = [1, 2, 3, 5, 10, 20, 50, 100, 200],
        gamma = 1.0,
        sigma = 1.0,
        epsilon = 0.0,
        broadening = 0.05,
        tmax = 80.0,
        nt = 3000,
        wmin = -5.0,
        wmax = 5.0,
        nw = 1201,
        max_hierarchy_size = max_hierarchy_size,
        filename = "identical_absorption_spectra.png"
    )

    make_identical_absorption_error_plot(
        nmin = 1,
        nmax = 500,
        npoints = 45,
        gamma = 1.0,
        sigma = 1.0,
        epsilon = 0.0,
        broadening = 0.05,
        tmax = 80.0,
        nt = 2500,
        wmin = -5.0,
        wmax = 5.0,
        nw = 801,
        max_hierarchy_size = max_hierarchy_size,
        filename = "identical_absorption_error_vs_n.png"
    )

    make_identical_absorption_gif(
        nmin = 1,
        nmax = 300,
        nframes = 40,
        fps = 8,
        gamma = 1.0,
        sigma = 1.0,
        epsilon = 0.0,
        broadening = 0.05,
        tmax = 80.0,
        nt = 2500,
        wmin = -5.0,
        wmax = 5.0,
        nw = 1001,
        max_hierarchy_size = max_hierarchy_size,
        filename = "identical_absorption_spectra.gif"
    )

    return nothing
end


function run_typed_examples(;
    K = 6,
    sigma = 1.0,
    gamma_min = 0.2,
    gamma_max = 5.0,
    amp_power = 0.35,
    amp_mod = 0.35,
    count_sampling = :balanced,
    gamma_sampling = :logrange,
    amplitude_sampling = :golden_modulated,
    count_weights = nothing,
    amplitude_weights = nothing,
    rng_seed = nothing,
    tag = "typed",
    max_hierarchy_size = DEFAULT_MAX_HIERARCHY_SIZE
)
    # For K = 6, the hierarchy size grows as prod(counts .+ 1).
    # These default n values are chosen so the exact hierarchy remains manageable.
    make_typed_distinct_absorption_plot(
        ns = [6, 12, 24, 36, 48, 60],
        K = K,
        sigma = sigma,
        gamma_min = gamma_min,
        gamma_max = gamma_max,
        amp_power = amp_power,
        amp_mod = amp_mod,
        count_sampling = count_sampling,
        gamma_sampling = gamma_sampling,
        amplitude_sampling = amplitude_sampling,
        count_weights = count_weights,
        amplitude_weights = amplitude_weights,
        rng_seed = rng_seed,
        epsilon = 0.0,
        broadening = 0.05,
        tmax = 80.0,
        nt = 3000,
        wmin = -5.0,
        wmax = 5.0,
        nw = 1201,
        max_hierarchy_size = max_hierarchy_size,
        filename = "$(tag)_absorption_spectra.png"
    )

    make_typed_distinct_absorption_error_plot(
        nmin = K,
        nmax = 60,
        npoints = 28,
        K = K,
        sigma = sigma,
        gamma_min = gamma_min,
        gamma_max = gamma_max,
        amp_power = amp_power,
        amp_mod = amp_mod,
        count_sampling = count_sampling,
        gamma_sampling = gamma_sampling,
        amplitude_sampling = amplitude_sampling,
        count_weights = count_weights,
        amplitude_weights = amplitude_weights,
        rng_seed = rng_seed,
        epsilon = 0.0,
        broadening = 0.05,
        tmax = 80.0,
        nt = 2500,
        wmin = -5.0,
        wmax = 5.0,
        nw = 801,
        max_hierarchy_size = max_hierarchy_size,
        filename = "$(tag)_absorption_error_vs_n.png"
    )

    make_typed_distinct_absorption_gif(
        nmin = K,
        nmax = 60,
        nframes = 30,
        fps = 8,
        K = K,
        sigma = sigma,
        gamma_min = gamma_min,
        gamma_max = gamma_max,
        amp_power = amp_power,
        amp_mod = amp_mod,
        count_sampling = count_sampling,
        gamma_sampling = gamma_sampling,
        amplitude_sampling = amplitude_sampling,
        count_weights = count_weights,
        amplitude_weights = amplitude_weights,
        rng_seed = rng_seed,
        epsilon = 0.0,
        broadening = 0.05,
        tmax = 80.0,
        nt = 2500,
        wmin = -5.0,
        wmax = 5.0,
        nw = 1001,
        max_hierarchy_size = max_hierarchy_size,
        filename = "$(tag)_absorption_spectra.gif"
    )

    return nothing
end


function run_unique_examples(;
    sigma = 1.0,
    gamma_min = 0.2,
    gamma_max = 5.0,
    amp_power = 0.35,
    amp_mod = 0.35,
    gamma_sampling = :logmidpoints,
    amplitude_sampling = :golden_modulated,
    amplitude_weights = nothing,
    rng_seed = nothing,
    tag = "unique",
    max_hierarchy_size = DEFAULT_MAX_HIERARCHY_SIZE
)
    # Exact unique-channel RTN-HE has hierarchy size 2^n.
    # Keeping n <= 20 is usually much safer for desktop-scale runs.
    make_unique_distinct_absorption_plot(
        ns = [2, 3, 4, 5, 6, 8, 10, 12],
        sigma = sigma,
        gamma_min = gamma_min,
        gamma_max = gamma_max,
        amp_power = amp_power,
        amp_mod = amp_mod,
        gamma_sampling = gamma_sampling,
        amplitude_sampling = amplitude_sampling,
        amplitude_weights = amplitude_weights,
        rng_seed = rng_seed,
        epsilon = 0.0,
        broadening = 0.05,
        tmax = 80.0,
        nt = 3000,
        wmin = -5.0,
        wmax = 5.0,
        nw = 1201,
        max_hierarchy_size = max_hierarchy_size,
        filename = "$(tag)_absorption_spectra.png"
    )

    make_unique_distinct_absorption_error_plot(
        nmin = 2,
        nmax = 20,
        npoints = 20,
        sigma = sigma,
        gamma_min = gamma_min,
        gamma_max = gamma_max,
        amp_power = amp_power,
        amp_mod = amp_mod,
        gamma_sampling = gamma_sampling,
        amplitude_sampling = amplitude_sampling,
        amplitude_weights = amplitude_weights,
        rng_seed = rng_seed,
        epsilon = 0.0,
        broadening = 0.05,
        tmax = 80.0,
        nt = 2500,
        wmin = -5.0,
        wmax = 5.0,
        nw = 801,
        max_hierarchy_size = max_hierarchy_size,
        filename = "$(tag)_absorption_error_vs_n.png"
    )

    make_unique_distinct_absorption_gif(
        nmin = 2,
        nmax = 18,
        nframes = 28,
        fps = 8,
        sigma = sigma,
        gamma_min = gamma_min,
        gamma_max = gamma_max,
        amp_power = amp_power,
        amp_mod = amp_mod,
        gamma_sampling = gamma_sampling,
        amplitude_sampling = amplitude_sampling,
        amplitude_weights = amplitude_weights,
        rng_seed = rng_seed,
        epsilon = 0.0,
        broadening = 0.05,
        tmax = 80.0,
        nt = 2500,
        wmin = -5.0,
        wmax = 5.0,
        nw = 1001,
        max_hierarchy_size = max_hierarchy_size,
        filename = "$(tag)_absorption_spectra.gif"
    )

    return nothing
end


function run_typed_random_examples()
    run_typed_examples(
        count_sampling = :random,
        gamma_sampling = :loguniform,
        amplitude_sampling = :lognormal,
        rng_seed = 1234,
        tag = "typed_random"
    )

    return nothing
end


function run_typed_weighted_examples()
    run_typed_examples(
        count_sampling = :weighted,
        count_weights = [1, 1, 2, 2, 4, 4],
        gamma_sampling = :loguniform,
        amplitude_sampling = :random_modulated,
        rng_seed = 1234,
        tag = "typed_weighted"
    )

    return nothing
end


function run_unique_random_examples()
    run_unique_examples(
        gamma_sampling = :loguniform,
        amplitude_sampling = :lognormal,
        rng_seed = 1234,
        tag = "unique_random"
    )

    return nothing
end


function run_all_examples()
    run_identical_examples()
    run_typed_examples()
    run_unique_examples()

    return nothing
end


function main(args = ARGS)
    mode = length(args) >= 1 ? lowercase(args[1]) : "identical"

    if mode == "identical"
        run_identical_examples()

    elseif mode == "typed"
        run_typed_examples()

    elseif mode in ("typed-random", "typed_random")
        run_typed_random_examples()

    elseif mode in ("typed-weighted", "typed_weighted")
        run_typed_weighted_examples()

    elseif mode == "unique"
        run_unique_examples()

    elseif mode in ("unique-random", "unique_random")
        run_unique_random_examples()

    elseif mode == "all"
        run_all_examples()

    else
        println("Unknown mode: $mode")
        println("Use one of: identical, typed, typed-random, typed-weighted, unique, unique-random, all")
    end

    return nothing
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
