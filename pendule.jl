include("imports.jl")

function degrees_to_radians(degrees)
    return degrees * (π / 180)
end

function Lagrangien(θ1_0, θ2_0, m1, m2, g, r1, r2, L1, L2)
    lagrange = m1*(L1*θ1_0)^2/2 + m2*((L1*θ1_0)^2 + (L2*θ2_0)^2 + 2*L1*L2*θ1_0*θ2_0*cos(θ1_0 - θ2_0)) + m1*g*L1*cos(θ1_0) + m2*g*(L1*cos(θ1_0) + L2*cos(θ2_0))
    return lagrange
end

function position(θ1_0, L1, θ2_0, L2)
    r1 = [L1*sin(θ1_0), L1*cos(θ1_0)]
    r2 = [r1[1] + L2*sin(θ2_0), r1[2] + L2*cos(θ2_0)]
    return r1, r2
end

# Parameters
L1 = 91.74  # Length of first pendulum (mm)
L2 = 69.33  # Length of second pendulum (mm)
θ1_0 = degrees_to_radians(181.5)  # Initial angle of first pendulum (rad)
θ2_0 = degrees_to_radians(183.1)  # Initial angle of second pendulum (rad)
g = 9.81    # Acceleration due to gravity (m/s^2)
m1 = 0      # Mass of first pendulum (to be determined)
m2 = 0      # Mass of second pendulum (to be determined)

# Initial positions
r1, r2 = position(θ1_0, L1, θ2_0, L2)
