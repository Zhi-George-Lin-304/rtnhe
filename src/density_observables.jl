# Utilities for reading physical information from density-matrix trajectories.

function density_populations(rhos::Array{ComplexF64, 3})
    d, d2, nt = size(rhos)
    @assert d == d2

    pops = zeros(Float64, nt, d)

    @inbounds for it in 1:nt
        for a in 1:d
            pops[it, a] = real(rhos[a, a, it])
        end
    end

    return pops
end

function density_coherence(rhos::Array{ComplexF64, 3}, i::Integer, j::Integer)
    d, d2, nt = size(rhos)
    @assert d == d2
    @assert 1 <= i <= d
    @assert 1 <= j <= d

    return [rhos[i, j, it] for it in 1:nt]
end

function density_expectation(rhos::Array{ComplexF64, 3}, O)
    d, d2, nt = size(rhos)
    @assert d == d2

    Oc = ComplexF64.(O)
    @assert size(Oc) == (d, d)

    vals = zeros(ComplexF64, nt)

    @inbounds for it in 1:nt
        rho = @view rhos[:, :, it]
        vals[it] = tr(Oc * rho)
    end

    return vals
end

function density_traces(rhos::Array{ComplexF64, 3})
    d, d2, nt = size(rhos)
    @assert d == d2

    vals = zeros(ComplexF64, nt)

    @inbounds for it in 1:nt
        vals[it] = tr(@view rhos[:, :, it])
    end

    return vals
end

function density_purities(rhos::Array{ComplexF64, 3})
    d, d2, nt = size(rhos)
    @assert d == d2

    vals = zeros(Float64, nt)

    @inbounds for it in 1:nt
        rho = @view rhos[:, :, it]
        vals[it] = real(tr(rho * rho))
    end

    return vals
end
