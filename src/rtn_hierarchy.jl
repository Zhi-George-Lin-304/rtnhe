using LinearAlgebra
using DifferentialEquations

# Store all parameters needed to propagate the grouped RTN hierarchy for optical coherence.
struct RTNGroupedCoherenceParams
    H0::Matrix{ComplexF64}
    P::Matrix{ComplexF64}
    gammas::Vector{Float64}
    deltas::Vector{Float64}
    counts::Vector{Int}
    dims::Vector{Int} # Hierarchy dimension for each group
    strides::Vector{Int}
    broadening::Float64
    d::Int # Hilbert-space dimension
    K::Int # Number of RTN groups
    m::Int # Total number of hierarchy components
end

# Build a parameter object for the grouped RTN optical-coherence hierarchy.
function make_rtn_grouped_coherence_params(
    H0,
    P,
    gammas,
    deltas,
    counts;
    broadening = 0.05
)
    K = length(counts)

    @assert length(gammas) == K
    @assert length(deltas) == K
    @assert all(counts .>= 1)
    @assert all(gammas .> 0)
    @assert broadening >= 0

    H0c = ComplexF64.(H0)
    Pc = ComplexF64.(P)

    d = size(H0c, 1)
    @assert size(H0c, 2) == d
    @assert size(Pc) == (d, d)

    # Group dimensions
    dims = counts .+ 1

    strides = ones(Int, K)
    for q in 2:K
        strides[q] = strides[q - 1] * dims[q - 1]
    end

    m = prod(dims)

    return RTNGroupedCoherenceParams(
        H0c,
        Pc,
        Float64.(gammas),
        Float64.(deltas),
        Int.(counts),
        Int.(dims),
        Int.(strides),
        Float64(broadening),
        d,
        K,
        m
    )
end

# Evaluate the right-hand side of the grouped RTN coherence hierarchy ODE.
function rtn_grouped_coherence_rhs!(du, u, p::RTNGroupedCoherenceParams, t)
    U = reshape(u, p.d, p.m)
    DU = reshape(du, p.d, p.m)

    @inbounds for idx in 1:p.m
        idx0 = idx - 1

        y = @view U[:, idx]
        dy = @view DU[:, idx]

        # Optical-coherence deterministic evolution: -i H0 y
        dy .= (-1im) .* (p.H0 * y)

        # Phenomenological homogeneous broadening / decay.
        if p.broadening != 0.0
            dy .-= p.broadening .* y
        end

        # Switching damping: -2 sum_q gamma_q r_q Y_r
        x = idx0
        damping = 0.0

        for q in 1:p.K
            rq = x % p.dims[q]
            x ÷= p.dims[q]
            damping += 2.0 * p.gammas[q] * rq
        end

        if damping != 0.0
            dy .-= damping .* y
        end

        # RTN hierarchy couplings.
        x = idx0

        for q in 1:p.K
            rq = x % p.dims[q]
            x ÷= p.dims[q]

            stride = p.strides[q]

            # Coupling to r_q - 1:
            # -i g_q P [r_q Y_{r-e_q}]
            if rq > 0
                ym = @view U[:, idx - stride]
                dy .+= (-1im * p.deltas[q] * rq) .* (p.P * ym)
            end

            # Coupling to r_q + 1:
            # -i g_q P [(m_q-r_q) Y_{r+e_q}]
            if rq < p.counts[q]
                yp = @view U[:, idx + stride]
                dy .+= (-1im * p.deltas[q] * (p.counts[q] - rq)) .* (p.P * yp)
            end
        end
    end

    return nothing
end

# Solve the grouped RTN hierarchy and return the time-domain optical response.
function solve_rtn_absorption_response(
    counts,
    gammas,
    deltas;
    H0,
    P,
    mu,
    broadening = 0.05,
    tmax = 80.0,
    nt = 3000,
    alg = Tsit5(),
    reltol = 1e-8,
    abstol = 1e-10
)
    params = make_rtn_grouped_coherence_params(
        H0,
        P,
        gammas,
        deltas,
        counts;
        broadening = broadening
    )

    muc = ComplexF64.(mu)
    @assert length(muc) == params.d

    U0 = zeros(ComplexF64, params.d, params.m)
    U0[:, 1] .= muc # Set the zeroth hierarchy component to the initial dipole vector.

    ts = collect(range(0.0, tmax, length = nt))

    prob = ODEProblem(
        rtn_grouped_coherence_rhs!,
        vec(U0),
        (0.0, tmax),
        params
    )

    sol = solve(
        prob,
        alg;
        saveat = ts,
        reltol = reltol,
        abstol = abstol
    )

    response = zeros(ComplexF64, length(sol.t))

    for i in eachindex(sol.u)
        U = reshape(sol.u[i], params.d, params.m)
        y0 = @view U[:, 1]
        response[i] = dot(muc, y0) 
    end

    return sol.t, response, sol, params
end