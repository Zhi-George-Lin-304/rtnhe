function gaussian_factor_grouped(t, counts, gammas, deltas)
    K = length(counts)
    exponent = 0.0

    @inbounds for q in 1:K
        gamma = gammas[q]
        delta = deltas[q]
        m = counts[q]

        exponent -= m * delta^2 * (
            t / (2.0 * gamma) -
            (1.0 - exp(-2.0 * gamma * t)) / (4.0 * gamma^2)
        )
    end

    return exp(exponent)
end

function gaussian_absorption_response(
    ts,
    counts,
    gammas,
    deltas;
    epsilon = 0.0,
    broadening = 0.05,
    dipole_strength = 1.0
)
    response = zeros(ComplexF64, length(ts))

    @inbounds for i in eachindex(ts)
        t = ts[i]
        response[i] = dipole_strength^2 *
                      exp((-1im * epsilon - broadening) * t) *
                      gaussian_factor_grouped(t, counts, gammas, deltas)
    end

    return response
end