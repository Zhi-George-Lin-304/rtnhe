
using Plots
using Printf
using Random

include("utils.jl")
include("simulation.jl")
include("channels.jl")


# ============================================================
# General helpers
# ============================================================

# Build a reproducible RNG if `rng_seed` is supplied.
_make_rng(rng_seed) = rng_seed === nothing ? Random.default_rng() : MersenneTwister(rng_seed)


# Estimate the RTN hierarchy size before solving.
# For grouped RTN-HE, the hierarchy size is prod(counts .+ 1).
function estimated_hierarchy_size(counts::AbstractVector{<:Integer})
    return prod(BigInt.(counts .+ 1))
end


# Skip cases whose hierarchy is too large to solve safely.
function _skip_if_too_large(counts, max_hierarchy_size; label = "")
    h = estimated_hierarchy_size(counts)
    hmax = BigInt(max_hierarchy_size)

    if h > hmax
        @warn "Skipping $label because estimated hierarchy size = $h exceeds max_hierarchy_size = $hmax."
        return true, h
    end

    return false, h
end


# Generate logarithmically spaced n values while enforcing a lower bound.
function _safe_logspace_ns(nmin, nmax, npoints; minimum_n = 1)
    ns = logspace_int(Int(nmin), Int(nmax), Int(npoints))
    ns = [max(Int(minimum_n), Int(n)) for n in ns]
    return unique(ns)
end


# Common call to solve_and_spectrum after the bath parameters are generated.
function _solve_checked_spectrum(
    counts,
    gammas,
    deltas;
    epsilon,
    broadening,
    tmax,
    nt,
    wmin,
    wmax,
    nw,
    max_hierarchy_size,
    label = ""
)
    should_skip, _ = _skip_if_too_large(
        counts,
        max_hierarchy_size;
        label = label
    )

    if should_skip
        return nothing
    end

    return solve_and_spectrum(
        counts,
        gammas,
        deltas;
        epsilon = epsilon,
        broadening = broadening,
        tmax = tmax,
        nt = nt,
        wmin = wmin,
        wmax = wmax,
        nw = nw
    )
end


# Plot one RTN curve and optionally the Gaussian benchmark.
function _add_absorption_curves!(
    p,
    result;
    label,
    add_gaussian = false,
    gaussian_label = "Gaussian",
    rtn_lw = 2,
    gaussian_lw = 4
)
    if add_gaussian
        plot!(
            p,
            result.omega,
            result.A_gaussian;
            lw = gaussian_lw,
            ls = :dash,
            label = gaussian_label
        )
    end

    plot!(
        p,
        result.omega,
        result.A_rtn;
        lw = rtn_lw,
        label = label
    )

    return p
end


# Standard error-vs-n plot.
function _make_error_plot(
    ns,
    errs;
    xlabel = "n",
    ylabel = "max_ω |A_RTN(ω) - A_Gauss(ω)|",
    title = "Absorption convergence to Gaussian",
    label = "spectral error",
    filename = "absorption_error_vs_n.png"
)
    @assert !isempty(ns) "No successful data points were generated."
    @assert length(ns) == length(errs)

    p = plot(
        ns,
        errs;
        marker = :circle,
        lw = 2,
        xscale = :log10,
        yscale = :log10,
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        label = label
    )

    hline!(p, [1e-2]; lw = 2, ls = :dash, label = "1e-2")
    hline!(p, [1e-3]; lw = 2, ls = :dot, label = "1e-3")

    savefig(p, filename)
    println("Saved $filename")

    return p
end


# ============================================================
# 1. Identical-channel RTN bath
# ============================================================

function make_identical_absorption_plot(;
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
    max_hierarchy_size = 2_000_000,
    filename = "identical_absorption_spectra.png"
)
    p = plot(
        xlabel = "frequency ω - ω₀",
        ylabel = "normalized absorption",
        title = "Absorption: identical RTN channels vs Gaussian",
        legend = :topright
    )

    println("Identical-channel absorption spectra:")

    added_gaussian = false

    for n in ns
        counts, gammas, deltas = identical_channels(n; gamma = gamma, sigma = sigma)

        result = _solve_checked_spectrum(
            counts,
            gammas,
            deltas;
            epsilon = epsilon,
            broadening = broadening,
            tmax = tmax,
            nt = nt,
            wmin = wmin,
            wmax = wmax,
            nw = nw,
            max_hierarchy_size = max_hierarchy_size,
            label = "identical n = $n"
        )

        result === nothing && continue

        _add_absorption_curves!(
            p,
            result;
            label = "RTN-HE n = $n",
            add_gaussian = !added_gaussian,
            gaussian_label = "Gaussian"
        )

        added_gaussian = true

        @printf(
            "n = %5d  hierarchy = %8d  n_eff = %8.2f  variance = %.4f  spectral error = %.4e\n",
            n,
            result.hierarchy_size,
            result.n_eff,
            total_noise_variance(result.counts, result.deltas),
            result.err
        )
    end

    savefig(p, filename)
    println("Saved $filename")

    return p
end


function make_identical_absorption_error_plot(;
    nmin = 1,
    nmax = 500,
    npoints = 50,
    gamma = 1.0,
    sigma = 1.0,
    epsilon = 0.0,
    broadening = 0.05,
    tmax = 80.0,
    nt = 2500,
    wmin = -5.0,
    wmax = 5.0,
    nw = 801,
    max_hierarchy_size = 2_000_000,
    filename = "identical_absorption_error_vs_n.png"
)
    ns_raw = _safe_logspace_ns(nmin, nmax, npoints; minimum_n = 1)
    ns = Int[]
    errs = Float64[]

    println("Identical-channel absorption error scan:")

    for n in ns_raw
        counts, gammas, deltas = identical_channels(n; gamma = gamma, sigma = sigma)

        result = _solve_checked_spectrum(
            counts,
            gammas,
            deltas;
            epsilon = epsilon,
            broadening = broadening,
            tmax = tmax,
            nt = nt,
            wmin = wmin,
            wmax = wmax,
            nw = nw,
            max_hierarchy_size = max_hierarchy_size,
            label = "identical n = $n"
        )

        result === nothing && continue

        push!(ns, n)
        push!(errs, result.err)

        @printf("n = %5d  spectral error = %.4e\n", n, result.err)
    end

    return _make_error_plot(
        ns,
        errs;
        title = "Identical RTN channels: convergence to Gaussian",
        filename = filename
    )
end


function make_identical_absorption_gif(;
    nmin = 1,
    nmax = 300,
    nframes = 45,
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
    max_hierarchy_size = 2_000_000,
    filename = "identical_absorption_spectra.gif"
)
    ns = _safe_logspace_ns(nmin, nmax, nframes; minimum_n = 1)

    anim = @animate for n in ns
        counts, gammas, deltas = identical_channels(n; gamma = gamma, sigma = sigma)

        result = _solve_checked_spectrum(
            counts,
            gammas,
            deltas;
            epsilon = epsilon,
            broadening = broadening,
            tmax = tmax,
            nt = nt,
            wmin = wmin,
            wmax = wmax,
            nw = nw,
            max_hierarchy_size = max_hierarchy_size,
            label = "identical n = $n"
        )

        if result === nothing
            plot(title = "Skipped n = $n: hierarchy too large", legend = false)
            continue
        end

        p = plot(
            result.omega,
            result.A_gaussian;
            lw = 4,
            ls = :dash,
            label = "Gaussian",
            xlabel = "frequency ω - ω₀",
            ylabel = "normalized absorption",
            title = "Identical RTN channels → Gaussian absorption",
            legend = :topright,
            ylim = (-0.05, 1.08)
        )

        plot!(
            p,
            result.omega,
            result.A_rtn;
            lw = 3,
            label = "RTN-HE"
        )

        txt = @sprintf(
            "n = %d\nhierarchy = %d\nn_eff = %.1f\nmax spectral diff = %.3e",
            n,
            result.hierarchy_size,
            result.n_eff,
            result.err
        )

        annotate!(p, wmin + 0.58 * (wmax - wmin), 0.28, text(txt, 11, :left))
        p
    end

    gif(anim, filename, fps = fps)
    println("Saved $filename")

    return filename
end


# Convenience wrapper: make all three identical-channel outputs.
# Use separate keyword bundles because PNG, error, and GIF functions have different keyword sets.
function make_all_identical_absorption_outputs(;
    plot_kwargs = NamedTuple(),
    error_kwargs = NamedTuple(),
    gif_kwargs = NamedTuple()
)
    p_png = make_identical_absorption_plot(; plot_kwargs...)
    p_err = make_identical_absorption_error_plot(; error_kwargs...)
    f_gif = make_identical_absorption_gif(; gif_kwargs...)
    return p_png, p_err, f_gif
end


# ============================================================
# 2. Typed distinct RTN groups
# ============================================================

function make_typed_distinct_absorption_plot(;
    ns = [6, 8, 10, 12, 16, 20, 24],
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
    epsilon = 0.0,
    broadening = 0.05,
    tmax = 80.0,
    nt = 3000,
    wmin = -5.0,
    wmax = 5.0,
    nw = 1201,
    max_hierarchy_size = 2_000_000,
    filename = "typed_distinct_absorption_spectra.png"
)
    rng = _make_rng(rng_seed)

    p = plot(
        xlabel = "frequency ω - ω₀",
        ylabel = "normalized absorption",
        title = "Absorption: typed RTN groups vs Gaussian",
        legend = :topright
    )

    println("Typed distinct-channel absorption spectra:")
    @printf(
        "K = %d, count_sampling = %s, gamma_sampling = %s, amplitude_sampling = %s\n",
        K,
        string(count_sampling),
        string(gamma_sampling),
        string(amplitude_sampling)
    )

    added_gaussian = false

    for n in ns
        counts, gammas, deltas = typed_distinct_groups(
            n;
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
            rng = rng
        )

        result = _solve_checked_spectrum(
            counts,
            gammas,
            deltas;
            epsilon = epsilon,
            broadening = broadening,
            tmax = tmax,
            nt = nt,
            wmin = wmin,
            wmax = wmax,
            nw = nw,
            max_hierarchy_size = max_hierarchy_size,
            label = "typed n = $n, K = $K"
        )

        result === nothing && continue

        _add_absorption_curves!(
            p,
            result;
            label = "RTN n = $n",
            add_gaussian = !added_gaussian,
            gaussian_label = "Gaussian"
        )

        added_gaussian = true

        @printf(
            "n = %5d  groups = %2d  hierarchy = %8d  n_eff = %8.2f  variance = %.4f  spectral error = %.4e\n",
            n,
            K,
            result.hierarchy_size,
            result.n_eff,
            total_noise_variance(result.counts, result.deltas),
            result.err
        )
    end

    savefig(p, filename)
    println("Saved $filename")

    return p
end


function make_typed_distinct_absorption_error_plot(;
    nmin = 6,
    nmax = 80,
    npoints = 30,
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
    epsilon = 0.0,
    broadening = 0.05,
    tmax = 80.0,
    nt = 2500,
    wmin = -5.0,
    wmax = 5.0,
    nw = 801,
    max_hierarchy_size = 2_000_000,
    filename = "typed_distinct_absorption_error_vs_n.png"
)
    rng = _make_rng(rng_seed)

    ns_raw = _safe_logspace_ns(nmin, nmax, npoints; minimum_n = K)
    ns = Int[]
    errs = Float64[]

    println("Typed distinct-channel absorption error scan:")
    @printf(
        "K = %d, count_sampling = %s, gamma_sampling = %s, amplitude_sampling = %s\n",
        K,
        string(count_sampling),
        string(gamma_sampling),
        string(amplitude_sampling)
    )

    for n in ns_raw
        counts, gammas, deltas = typed_distinct_groups(
            n;
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
            rng = rng
        )

        result = _solve_checked_spectrum(
            counts,
            gammas,
            deltas;
            epsilon = epsilon,
            broadening = broadening,
            tmax = tmax,
            nt = nt,
            wmin = wmin,
            wmax = wmax,
            nw = nw,
            max_hierarchy_size = max_hierarchy_size,
            label = "typed n = $n, K = $K"
        )

        result === nothing && continue

        push!(ns, n)
        push!(errs, result.err)

        @printf(
            "n = %5d  groups = %2d  hierarchy = %8d  spectral error = %.4e\n",
            n,
            K,
            result.hierarchy_size,
            result.err
        )
    end

    return _make_error_plot(
        ns,
        errs;
        title = "Typed distinct RTN groups: convergence to Gaussian",
        filename = filename
    )
end


function make_typed_distinct_absorption_gif(;
    nmin = 6,
    nmax = 80,
    nframes = 35,
    fps = 8,
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
    epsilon = 0.0,
    broadening = 0.05,
    tmax = 80.0,
    nt = 2500,
    wmin = -5.0,
    wmax = 5.0,
    nw = 1001,
    max_hierarchy_size = 2_000_000,
    filename = "typed_distinct_absorption_spectra.gif"
)
    rng = _make_rng(rng_seed)

    ns = _safe_logspace_ns(nmin, nmax, nframes; minimum_n = K)

    anim = @animate for n in ns
        counts, gammas, deltas = typed_distinct_groups(
            n;
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
            rng = rng
        )

        result = _solve_checked_spectrum(
            counts,
            gammas,
            deltas;
            epsilon = epsilon,
            broadening = broadening,
            tmax = tmax,
            nt = nt,
            wmin = wmin,
            wmax = wmax,
            nw = nw,
            max_hierarchy_size = max_hierarchy_size,
            label = "typed n = $n, K = $K"
        )

        if result === nothing
            plot(title = "Skipped n = $n: hierarchy too large", legend = false)
            continue
        end

        p = plot(
            result.omega,
            result.A_gaussian;
            lw = 4,
            ls = :dash,
            label = "Gaussian with same covariance",
            xlabel = "frequency ω - ω₀",
            ylabel = "normalized absorption",
            title = "Typed distinct RTN groups → Gaussian absorption",
            legend = :topright,
            ylim = (-0.05, 1.08)
        )

        plot!(
            p,
            result.omega,
            result.A_rtn;
            lw = 3,
            label = "direct RTN-HE"
        )

        txt = @sprintf(
            "n = %d, K = %d\n%s / %s / %s\nhierarchy = %d\nn_eff = %.1f\nmax diff = %.3e",
            n,
            K,
            string(count_sampling),
            string(gamma_sampling),
            string(amplitude_sampling),
            result.hierarchy_size,
            result.n_eff,
            result.err
        )

        annotate!(p, wmin + 0.50 * (wmax - wmin), 0.30, text(txt, 10, :left))
        p
    end

    gif(anim, filename, fps = fps)
    println("Saved $filename")

    return filename
end


# Convenience wrapper: make all three typed-group outputs.
# Use separate keyword bundles because PNG, error, and GIF functions have different keyword sets.
function make_all_typed_distinct_absorption_outputs(;
    plot_kwargs = NamedTuple(),
    error_kwargs = NamedTuple(),
    gif_kwargs = NamedTuple()
)
    p_png = make_typed_distinct_absorption_plot(; plot_kwargs...)
    p_err = make_typed_distinct_absorption_error_plot(; error_kwargs...)
    f_gif = make_typed_distinct_absorption_gif(; gif_kwargs...)
    return p_png, p_err, f_gif
end


# ============================================================
# 3. Unique distinct RTN channels
# ============================================================

function make_unique_distinct_absorption_plot(;
    ns = [2, 3, 4, 5, 6, 8, 10, 12],
    sigma = 1.0,
    gamma_min = 0.2,
    gamma_max = 5.0,
    amp_power = 0.35,
    amp_mod = 0.35,
    gamma_sampling = :logmidpoints,
    amplitude_sampling = :golden_modulated,
    amplitude_weights = nothing,
    rng_seed = nothing,
    epsilon = 0.0,
    broadening = 0.05,
    tmax = 80.0,
    nt = 3000,
    wmin = -5.0,
    wmax = 5.0,
    nw = 1201,
    max_hierarchy_size = 2_000_000,
    filename = "unique_distinct_absorption_spectra.png"
)
    rng = _make_rng(rng_seed)

    p = plot(
        xlabel = "frequency ω - ω₀",
        ylabel = "normalized absorption",
        title = "Absorption: unique distinct RTN channels",
        legend = :topright
    )

    println("Unique distinct-channel absorption spectra:")
    println("Warning: hierarchy size is 2^n for this exact unique-channel case.")
    @printf(
        "gamma_sampling = %s, amplitude_sampling = %s\n",
        string(gamma_sampling),
        string(amplitude_sampling)
    )

    added_gaussian = false

    for n in ns
        counts, gammas, deltas = unique_distinct_channels(
            n;
            sigma = sigma,
            gamma_min = gamma_min,
            gamma_max = gamma_max,
            amp_power = amp_power,
            amp_mod = amp_mod,
            gamma_sampling = gamma_sampling,
            amplitude_sampling = amplitude_sampling,
            amplitude_weights = amplitude_weights,
            rng = rng
        )

        result = _solve_checked_spectrum(
            counts,
            gammas,
            deltas;
            epsilon = epsilon,
            broadening = broadening,
            tmax = tmax,
            nt = nt,
            wmin = wmin,
            wmax = wmax,
            nw = nw,
            max_hierarchy_size = max_hierarchy_size,
            label = "unique n = $n"
        )

        result === nothing && continue

        _add_absorption_curves!(
            p,
            result;
            label = "RTN n = $n",
            add_gaussian = !added_gaussian,
            gaussian_label = "Gaussian"
        )

        added_gaussian = true

        @printf(
            "n = %4d  hierarchy = %8d  n_eff = %8.2f  variance = %.4f  spectral error = %.4e\n",
            n,
            result.hierarchy_size,
            result.n_eff,
            total_noise_variance(result.counts, result.deltas),
            result.err
        )
    end

    savefig(p, filename)
    println("Saved $filename")

    return p
end


function make_unique_distinct_absorption_error_plot(;
    nmin = 2,
    nmax = 20,
    npoints = 20,
    sigma = 1.0,
    gamma_min = 0.2,
    gamma_max = 5.0,
    amp_power = 0.35,
    amp_mod = 0.35,
    gamma_sampling = :logmidpoints,
    amplitude_sampling = :golden_modulated,
    amplitude_weights = nothing,
    rng_seed = nothing,
    epsilon = 0.0,
    broadening = 0.05,
    tmax = 80.0,
    nt = 2500,
    wmin = -5.0,
    wmax = 5.0,
    nw = 801,
    max_hierarchy_size = 2_000_000,
    filename = "unique_distinct_absorption_error_vs_n.png"
)
    rng = _make_rng(rng_seed)

    ns_raw = _safe_logspace_ns(nmin, nmax, npoints; minimum_n = 1)
    ns = Int[]
    errs = Float64[]

    println("Unique distinct-channel absorption error scan:")
    println("Warning: hierarchy size is 2^n for this exact unique-channel case.")
    @printf(
        "gamma_sampling = %s, amplitude_sampling = %s\n",
        string(gamma_sampling),
        string(amplitude_sampling)
    )

    for n in ns_raw
        counts, gammas, deltas = unique_distinct_channels(
            n;
            sigma = sigma,
            gamma_min = gamma_min,
            gamma_max = gamma_max,
            amp_power = amp_power,
            amp_mod = amp_mod,
            gamma_sampling = gamma_sampling,
            amplitude_sampling = amplitude_sampling,
            amplitude_weights = amplitude_weights,
            rng = rng
        )

        result = _solve_checked_spectrum(
            counts,
            gammas,
            deltas;
            epsilon = epsilon,
            broadening = broadening,
            tmax = tmax,
            nt = nt,
            wmin = wmin,
            wmax = wmax,
            nw = nw,
            max_hierarchy_size = max_hierarchy_size,
            label = "unique n = $n"
        )

        result === nothing && continue

        push!(ns, n)
        push!(errs, result.err)

        @printf(
            "n = %4d  hierarchy = %8d  spectral error = %.4e\n",
            n,
            result.hierarchy_size,
            result.err
        )
    end

    return _make_error_plot(
        ns,
        errs;
        title = "Unique distinct RTN channels: convergence to Gaussian",
        filename = filename
    )
end


function make_unique_distinct_absorption_gif(;
    nmin = 2,
    nmax = 18,
    nframes = 30,
    fps = 8,
    sigma = 1.0,
    gamma_min = 0.2,
    gamma_max = 5.0,
    amp_power = 0.35,
    amp_mod = 0.35,
    gamma_sampling = :logmidpoints,
    amplitude_sampling = :golden_modulated,
    amplitude_weights = nothing,
    rng_seed = nothing,
    epsilon = 0.0,
    broadening = 0.05,
    tmax = 80.0,
    nt = 2500,
    wmin = -5.0,
    wmax = 5.0,
    nw = 1001,
    max_hierarchy_size = 2_000_000,
    filename = "unique_distinct_absorption_spectra.gif"
)
    rng = _make_rng(rng_seed)

    ns = _safe_logspace_ns(nmin, nmax, nframes; minimum_n = 1)

    anim = @animate for n in ns
        counts, gammas, deltas = unique_distinct_channels(
            n;
            sigma = sigma,
            gamma_min = gamma_min,
            gamma_max = gamma_max,
            amp_power = amp_power,
            amp_mod = amp_mod,
            gamma_sampling = gamma_sampling,
            amplitude_sampling = amplitude_sampling,
            amplitude_weights = amplitude_weights,
            rng = rng
        )

        result = _solve_checked_spectrum(
            counts,
            gammas,
            deltas;
            epsilon = epsilon,
            broadening = broadening,
            tmax = tmax,
            nt = nt,
            wmin = wmin,
            wmax = wmax,
            nw = nw,
            max_hierarchy_size = max_hierarchy_size,
            label = "unique n = $n"
        )

        if result === nothing
            plot(title = "Skipped n = $n: hierarchy too large", legend = false)
            continue
        end

        p = plot(
            result.omega,
            result.A_gaussian;
            lw = 4,
            ls = :dash,
            label = "Gaussian with same covariance",
            xlabel = "frequency ω - ω₀",
            ylabel = "normalized absorption",
            title = "Unique distinct RTN channels → Gaussian absorption",
            legend = :topright,
            ylim = (-0.05, 1.08)
        )

        plot!(
            p,
            result.omega,
            result.A_rtn;
            lw = 3,
            label = "direct RTN-HE"
        )

        txt = @sprintf(
            "n = %d\n%s / %s\nhierarchy = %d\nn_eff = %.1f\nmax diff = %.3e",
            n,
            string(gamma_sampling),
            string(amplitude_sampling),
            result.hierarchy_size,
            result.n_eff,
            result.err
        )

        annotate!(p, wmin + 0.52 * (wmax - wmin), 0.30, text(txt, 10, :left))
        p
    end

    gif(anim, filename, fps = fps)
    println("Saved $filename")

    return filename
end


# Convenience wrapper: make all three unique-channel outputs.
# Use separate keyword bundles because PNG, error, and GIF functions have different keyword sets.
function make_all_unique_distinct_absorption_outputs(;
    plot_kwargs = NamedTuple(),
    error_kwargs = NamedTuple(),
    gif_kwargs = NamedTuple()
)
    p_png = make_unique_distinct_absorption_plot(; plot_kwargs...)
    p_err = make_unique_distinct_absorption_error_plot(; error_kwargs...)
    f_gif = make_unique_distinct_absorption_gif(; gif_kwargs...)
    return p_png, p_err, f_gif
end


# ============================================================
# Example calls
# ============================================================

# Identical channels:
#
# make_identical_absorption_plot()
# make_identical_absorption_error_plot()
# make_identical_absorption_gif()
#
# Or:
#
# make_all_identical_absorption_outputs()


# Typed distinct groups with deterministic original-like sampling:
#
# make_typed_distinct_absorption_plot(
#     K = 6,
#     count_sampling = :balanced,
#     gamma_sampling = :logrange,
#     amplitude_sampling = :golden_modulated
# )
#
# make_typed_distinct_absorption_error_plot(
#     K = 6,
#     count_sampling = :balanced,
#     gamma_sampling = :logrange,
#     amplitude_sampling = :golden_modulated
# )
#
# make_typed_distinct_absorption_gif(
#     K = 6,
#     count_sampling = :balanced,
#     gamma_sampling = :logrange,
#     amplitude_sampling = :golden_modulated
# )


# Typed distinct groups with random sampling:
#
# make_all_typed_distinct_absorption_outputs(
#     plot_kwargs = (
#         K = 6,
#         count_sampling = :random,
#         gamma_sampling = :loguniform,
#         amplitude_sampling = :lognormal,
#         rng_seed = 1234,
#     ),
#     error_kwargs = (
#         K = 6,
#         count_sampling = :random,
#         gamma_sampling = :loguniform,
#         amplitude_sampling = :lognormal,
#         rng_seed = 1234,
#     ),
#     gif_kwargs = (
#         K = 6,
#         count_sampling = :random,
#         gamma_sampling = :loguniform,
#         amplitude_sampling = :lognormal,
#         rng_seed = 1234,
#     )
# )


# Typed distinct groups with weighted group sizes:
#
# make_all_typed_distinct_absorption_outputs(
#     plot_kwargs = (
#         K = 6,
#         count_sampling = :weighted,
#         count_weights = [1, 1, 2, 2, 4, 4],
#         gamma_sampling = :loguniform,
#         amplitude_sampling = :random_modulated,
#         rng_seed = 1234,
#     ),
#     error_kwargs = (
#         K = 6,
#         count_sampling = :weighted,
#         count_weights = [1, 1, 2, 2, 4, 4],
#         gamma_sampling = :loguniform,
#         amplitude_sampling = :random_modulated,
#         rng_seed = 1234,
#     ),
#     gif_kwargs = (
#         K = 6,
#         count_sampling = :weighted,
#         count_weights = [1, 1, 2, 2, 4, 4],
#         gamma_sampling = :loguniform,
#         amplitude_sampling = :random_modulated,
#         rng_seed = 1234,
#     )
# )


# Unique distinct channels:
#
# make_unique_distinct_absorption_plot()
# make_unique_distinct_absorption_error_plot()
# make_unique_distinct_absorption_gif()
#
# Or:
#
# make_all_unique_distinct_absorption_outputs(
#     plot_kwargs = (
#         gamma_sampling = :loguniform,
#         amplitude_sampling = :lognormal,
#         rng_seed = 1234,
#     ),
#     error_kwargs = (
#         gamma_sampling = :loguniform,
#         amplitude_sampling = :lognormal,
#         rng_seed = 1234,
#     ),
#     gif_kwargs = (
#         gamma_sampling = :loguniform,
#         amplitude_sampling = :lognormal,
#         rng_seed = 1234,
#     )
# )
