# Construct the minimal one-level optical transition model
# used for RTN absorption-spectrum calculations.
function default_single_transition_model(; epsilon = 0.0)
    H0 = reshape(ComplexF64[epsilon], 1, 1)
    P = reshape(ComplexF64[1.0], 1, 1)
    mu = ComplexF64[1.0]
    return H0, P, mu
end