include("imports.jl")

md"""
# Double pendulum simulation

This is a simulation of a double pendulum using Julia.

The objectif is too recreate a double pendulum simulation based on videos and to predict his mouvement during 2 seconds past the video end.
"""

md"""
## Theory
A double pendulum is also know as Chaotic pendulum because of its sensibility to initial conditions.

The particularity of this system is his behavior. It can vary from simple periodic motion to a more complex and chaotic one with a small change of initial angles or velocities.

"""

md"""
## length of the pendulums (via gimps)
pendulum one = 91.74 mm
pendulum two = 69.33 mm
"""

md"""
## Parameters
- length of the first pendulum: L1 = 91.74 mm
- length of the second pendulum: L2 =  69.33 mm
- angles at t=0: θ1 = 181.5° = 3.159, θ2 = 183.1° = 3.194 (θ = 0° is the vertical down position)
"""

md"""
## Energy loss
In a real system, energy is lost due to friction and air resistance. This can be modeled by adding a damping term to the equations of motion.
We will have something like this:
$\begin{pmatrix}m_1\\m_2\\\alpha_1\\\alpha_2\end{pmatrix}$
where α1 and α2 are the angular accelerations of the first and second pendulums, respectively.
They will represent the energy loss in the system.

This energy loss is negligible for the first 1 second of the video, wich will help to find the mass of the two pendulums.
From the second 1, the energy loss will be taken into account to predict the motion of the pendulum after the video end.

The approch will be to find the mass of the two pendulums that minimize the difference between the simulated and observed positions over the first second.
And to find the angular accelerations, one of the approch will be to calculate the difference of energy between two frames and to divide it by the moment of inertia of each pendulum.
Or to test different values of angular accelerations to see which one fit the best with the observed data.
"""

function degrees_to_radians(degrees)
    return degrees * (π / 180)
end

function Lagrangien(θ1_0, θ2_0, m1, m2, g, r1, r2)
    # TODO: Implement Lagrangian mechanics
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
