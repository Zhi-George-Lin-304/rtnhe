# Sample a categorical index from probabilities `probs`.
function _sample_categorical(rng, probs::AbstractVector{<:Real})
    u = rand(rng)
    s = 0.0

    @inbounds for i in eachindex(probs)
        s += probs[i]
        if u <= s
            return i
        end
    end

    return length(probs)
end


# Normalize RTN amplitudes so that the total variance satisfies
#     sum_q counts[q] * deltas[q]^2 = sigma^2.
function normalize_rtn_deltas(
    counts::AbstractVector{<:Integer},
    weights::AbstractVector{<:Real};
    sigma = 1.0
)
    @assert length(counts) == length(weights)
    @assert sigma >= 0
    @assert all(counts .>= 1)
    @assert all(weights .>= 0)

    denom = sqrt(sum(counts .* weights.^2))
    @assert denom > 0 "At least one amplitude weight must be positive."

    deltas = sigma .* weights ./ denom

    return deltas
end


# Sample how many RTN channels belong to each group.
function sample_rtn_counts(
    n::Integer,
    K::Integer;
    method = :balanced,
    rng = Random.default_rng(),
    count_weights = nothing
)
    @assert K >= 1
    @assert n >= K "For grouped channels, choose n >= K."

    if method == :balanced
        counts = fill(div(n, K), K)
        remainder = n - sum(counts)

        for q in 1:remainder
            counts[q] += 1
        end

    elseif method == :random
        counts = ones(Int, K)

        for _ in 1:(n - K)
            q = rand(rng, 1:K)
            counts[q] += 1
        end

    elseif method == :weighted
        @assert count_weights !== nothing "Provide `count_weights` when method = :weighted."
        @assert length(count_weights) == K
        @assert all(count_weights .>= 0)
        @assert sum(count_weights) > 0

        probs = Float64.(count_weights) ./ sum(count_weights)
        counts = ones(Int, K)

        for _ in 1:(n - K)
            q = _sample_categorical(rng, probs)
            counts[q] += 1
        end

    else
        error("Unknown count sampling method: $method")
    end

    return counts
end


# Sample RTN switching rates gamma_q.
function sample_rtn_gammas(
    K::Integer;
    method = :logrange,
    gamma = 1.0,
    gamma_min = 0.2,
    gamma_max = 5.0,
    rng = Random.default_rng(),
    sort_values = true
)
    @assert K >= 1
    @assert gamma > 0
    @assert gamma_min > 0
    @assert gamma_max > 0
    @assert gamma_max >= gamma_min

    if method == :constant
        gammas = fill(Float64(gamma), K)

    elseif method == :logrange
        if K == 1
            gammas = [sqrt(gamma_min * gamma_max)]
        else
            gammas = collect(exp.(range(log(gamma_min), log(gamma_max), length = K)))
        end

    elseif method == :logmidpoints
        u = (collect(1:K) .- 0.5) ./ K
        gammas = exp.(log(gamma_min) .+ u .* (log(gamma_max) - log(gamma_min)))

    elseif method == :linrange
        if K == 1
            gammas = [(gamma_min + gamma_max) / 2]
        else
            gammas = collect(range(gamma_min, gamma_max, length = K))
        end

    elseif method == :linmidpoints
        u = (collect(1:K) .- 0.5) ./ K
        gammas = gamma_min .+ u .* (gamma_max - gamma_min)

    elseif method == :loguniform
        gammas = exp.(log(gamma_min) .+ rand(rng, K) .* (log(gamma_max) - log(gamma_min)))

        if sort_values
            sort!(gammas)
        end

    elseif method == :uniform
        gammas = gamma_min .+ rand(rng, K) .* (gamma_max - gamma_min)

        if sort_values
            sort!(gammas)
        end

    else
        error("Unknown gamma sampling method: $method")
    end

    return gammas
end


# Sample unnormalized RTN amplitude weights before variance normalization.
function sample_rtn_amplitude_weights(
    gammas::AbstractVector{<:Real};
    method = :golden_modulated,
    gamma_min = minimum(gammas),
    gamma_max = maximum(gammas),
    amp_power = 0.35,
    amp_mod = 0.35,
    rng = Random.default_rng(),
    amplitude_weights = nothing
)
    K = length(gammas)

    @assert K >= 1
    @assert all(gammas .> 0)
    @assert gamma_min > 0
    @assert gamma_max > 0
    @assert gamma_max >= gamma_min

    gamma_center = sqrt(gamma_min * gamma_max)
    gamma_factor = (gammas ./ gamma_center).^amp_power

    if method == :equal
        weights = ones(Float64, K)

    elseif method == :gamma_power
        weights = gamma_factor

    elseif method == :golden_modulated
        idx = collect(1:K)
        golden = (sqrt(5.0) - 1.0) / 2.0

        weights =
            gamma_factor .*
            exp.(amp_mod .* sin.(2.0 * pi .* golden .* idx))

    elseif method == :random_modulated
        weights =
            gamma_factor .*
            exp.(amp_mod .* (2.0 .* rand(rng, K) .- 1.0))

    elseif method == :lognormal
        weights =
            gamma_factor .*
            exp.(amp_mod .* randn(rng, K))

    elseif method == :given
        @assert amplitude_weights !== nothing "Provide `amplitude_weights` when method = :given."
        @assert length(amplitude_weights) == K
        weights = Float64.(amplitude_weights)

    else
        error("Unknown amplitude sampling method: $method")
    end

    @assert all(isfinite.(weights))
    @assert all(weights .>= 0)
    @assert sum(weights) > 0

    return weights
end


# Generate an identical-channel RTN bath with identical switching rates
# and couplings scaled to maintain fixed total variance.
function identical_channels(
    n::Integer;
    gamma = 1.0,
    sigma = 1.0
)
    @assert n >= 1
    @assert gamma > 0
    @assert sigma >= 0

    counts = [Int(n)]
    gammas = [Float64(gamma)]
    deltas = [Float64(sigma) / sqrt(n)]

    return counts, gammas, deltas
end


# Generate a heterogeneous RTN bath composed of K groups of fluctuators
# with selectable count, switching-rate, and amplitude sampling methods.
function typed_distinct_groups(
    n::Integer;
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
    rng = Random.default_rng()
)
    @assert K >= 1
    @assert n >= K "For this helper choose n >= K."
    @assert sigma >= 0
    @assert gamma_min > 0
    @assert gamma_max > 0
    @assert gamma_max >= gamma_min

    counts = sample_rtn_counts(
        n,
        K;
        method = count_sampling,
        rng = rng,
        count_weights = count_weights
    )

    gammas = sample_rtn_gammas(
        K;
        method = gamma_sampling,
        gamma_min = gamma_min,
        gamma_max = gamma_max,
        rng = rng
    )

    weights = sample_rtn_amplitude_weights(
        gammas;
        method = amplitude_sampling,
        gamma_min = gamma_min,
        gamma_max = gamma_max,
        amp_power = amp_power,
        amp_mod = amp_mod,
        rng = rng,
        amplitude_weights = amplitude_weights
    )

    deltas = normalize_rtn_deltas(
        counts,
        weights;
        sigma = sigma
    )

    return counts, gammas, deltas
end


# Generate a fully heterogeneous RTN bath in which every channel has its own switching rate
# and coupling strength, with selectable switching-rate and amplitude sampling methods.
function unique_distinct_channels(
    n::Integer;
    sigma = 1.0,
    gamma_min = 0.2,
    gamma_max = 5.0,
    amp_power = 0.35,
    amp_mod = 0.35,
    gamma_sampling = :logmidpoints,
    amplitude_sampling = :golden_modulated,
    amplitude_weights = nothing,
    rng = Random.default_rng()
)
    @assert n >= 1
    @assert sigma >= 0
    @assert gamma_min > 0
    @assert gamma_max > 0
    @assert gamma_max >= gamma_min

    counts = ones(Int, n)

    gammas = sample_rtn_gammas(
        n;
        method = gamma_sampling,
        gamma_min = gamma_min,
        gamma_max = gamma_max,
        rng = rng
    )

    weights = sample_rtn_amplitude_weights(
        gammas;
        method = amplitude_sampling,
        gamma_min = gamma_min,
        gamma_max = gamma_max,
        amp_power = amp_power,
        amp_mod = amp_mod,
        rng = rng,
        amplitude_weights = amplitude_weights
    )

    deltas = normalize_rtn_deltas(
        counts,
        weights;
        sigma = sigma
    )

    return counts, gammas, deltas
end