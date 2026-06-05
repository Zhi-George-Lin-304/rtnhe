# Return approximately logarithmically spaced integer values between `nmin` and
# `nmax`, using `m` sample points before rounding and duplicate removal.
function logspace_int(nmin::Int, nmax::Int, m::Int)
    @assert nmin >= 1
    @assert nmax >= nmin
    vals = round.(Int, exp.(range(log(Float64(nmin)), log(Float64(nmax)), length = m)))
    return unique(vals)
end

# Construct trapezoidal-rule quadrature weights for a one-dimensional time grid
# `ts`, including the case of a nonuniform grid.
function trapezoid_weights(ts::AbstractVector{<:Real})
    n = length(ts)
    @assert n >= 2

    w = zeros(Float64, n)
    w[1] = 0.5 * (ts[2] - ts[1])
    w[end] = 0.5 * (ts[end] - ts[end - 1])

    for i in 2:(n - 1)
        w[i] = 0.5 * (ts[i + 1] - ts[i - 1])
    end

    return w
end

# Normalize a real-valued vector in place by dividing by its maximum absolute
# value, while safely handling the all-zero vector.
function safe_normalize!(A::Vector{Float64})
    m = maximum(abs.(A))
    if m > 0
        A ./= m
    end
    return A
end

# Compute the n_eff of RTN channels.
function effective_channel_number(counts::Vector{Int}, gs::Vector{Float64})
    s2 = sum(counts .* gs.^2)
    s4 = sum(counts .* gs.^4)
    return s2^2 / s4
end

# Compute the total noise variance of RTN channels.
function total_noise_variance(counts::Vector{Int}, gs::Vector{Float64})
    return sum(counts .* gs.^2)
end