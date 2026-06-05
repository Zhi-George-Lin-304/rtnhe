include("models.jl")
include("rtn_hierarchy.jl")
include("gaussian_benchmark.jl")
include("spectrum.jl")

function solve_and_spectrum(
    counts,
    gammas,
    deltas;
    epsilon = 0.0,
    broadening = 0.05,
    tmax = 80.0,
    nt = 3000,
    wmin = -6.0,
    wmax = 6.0,
    nw = 1201,
    normalize = true,
    alg = Tsit5()
)
    H0, P, mu = default_single_transition_model(epsilon = epsilon)

    ts, response_rtn, sol, params = solve_rtn_absorption_response(
        counts,
        gammas,
        deltas;
        H0 = H0,
        P = P,
        mu = mu,
        broadening = broadening,
        tmax = tmax,
        nt = nt,
        alg = alg
    )

    response_g = gaussian_absorption_response(
        ts,
        counts,
        gammas,
        deltas;
        epsilon = epsilon,
        broadening = broadening,
        dipole_strength = abs(mu[1])
    )

    ws, A_rtn = absorption_spectrum(
        ts,
        response_rtn;
        wmin = wmin,
        wmax = wmax,
        nw = nw,
        normalize = normalize
    )

    _, A_g = absorption_spectrum(
        ts,
        response_g;
        wmin = wmin,
        wmax = wmax,
        nw = nw,
        normalize = normalize
    )

    err = maximum(abs.(A_rtn .- A_g))

    return (
        t = ts,
        response_rtn = response_rtn,
        response_gaussian = response_g,
        omega = ws,
        A_rtn = A_rtn,
        A_gaussian = A_g,
        err = err,
        n_eff = effective_channel_number(Int.(counts), Float64.(deltas)),
        hierarchy_size = params.m,
        counts = Int.(counts),
        gammas = Float64.(gammas),
        deltas = Float64.(deltas),
        sol = sol,
        params = params
    )
end