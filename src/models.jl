# Construct the minimal one-level optical transition model
# used for RTN absorption-spectrum calculations.
function default_single_transition_model(; epsilon = 0.0)
    H0 = reshape(ComplexF64[epsilon], 1, 1)
    P = reshape(ComplexF64[1.0], 1, 1)
    mu = ComplexF64[1.0]
    return H0, P, mu
end

function default_two_level_density_model(;
    epsilon = 1.0,
    J = 0.0,
    noise_on = :excited,
    initial = :plus
)
    H0 = ComplexF64[
        0.0      J;
        J  epsilon
    ]

    if noise_on == :excited
        P = ComplexF64[
            0.0  0.0;
            0.0  1.0
        ]
    elseif noise_on == :ground
        P = ComplexF64[
            1.0  0.0;
            0.0  0.0
        ]
    elseif noise_on == :sigma_z
        P = ComplexF64[
            0.5   0.0;
            0.0  -0.5
        ]
    else
        error("Unknown noise_on value: $noise_on. Use :excited, :ground, or :sigma_z.")
    end

    if initial == :ground
        psi = ComplexF64[1.0, 0.0]
    elseif initial == :excited
        psi = ComplexF64[0.0, 1.0]
    elseif initial == :plus
        psi = ComplexF64[1.0, 1.0] ./ sqrt(2.0)
    elseif initial == :minus
        psi = ComplexF64[1.0, -1.0] ./ sqrt(2.0)
    else
        error("Unknown initial value: $initial. Use :ground, :excited, :plus, or :minus.")
    end

    rho0 = psi * psi'

    return H0, P, rho0
end