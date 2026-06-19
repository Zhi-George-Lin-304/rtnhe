using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include(joinpath(@__DIR__, "..", "src", "RTNHE.jl"))
using .RTNHE
using Plots

gr()

counts, gammas, deltas = identical_channels(20; gamma = 1.0, sigma = 1.0)

H0, P, rho0 = default_two_level_density_model(
    epsilon = 1.0,
    J = 0.2,
    initial = :plus
)

ts, rhos, sol, params = solve_rtn_density_matrices(
    counts,
    gammas,
    deltas;
    H0 = H0,
    P = P,
    rho0 = rho0,
    tmax = 30.0,
    nt = 1000
)

pops = density_populations(rhos)
coh12 = density_coherence(rhos, 1, 2)

p = plot(
    ts,
    pops[:, 1],
    label = "rho11",
    xlabel = "time",
    ylabel = "density-matrix element"
)

plot!(p, ts, pops[:, 2], label = "rho22")
plot!(p, ts, abs.(coh12), label = "|rho12|")

savefig(p, "density_matrix_demo.png")