include("imports.jl")

# =============================================================================
# Chargement des données Tracker (fichier .trk)
# =============================================================================

function load_tracker_positions(filepath::String)
    """Charge les positions des masses depuis le fichier .trk"""
    xdoc = parse_file(filepath)
    xroot = root(xdoc)
    
    origin_x = 0.0
    origin_y = 0.0
    scale = 1.0
    
    for prop in child_elements(xroot)
        if attribute(prop, "name") == "coords"
            for obj in child_elements(prop)
                for subprop in child_elements(obj)
                    if attribute(subprop, "name") == "framedata"
                        for arr in child_elements(subprop)
                            for framedata in child_elements(arr)
                                for fprop in child_elements(framedata)
                                    name = attribute(fprop, "name")
                                    if name == "xorigin"
                                        origin_x = parse(Float64, content(fprop))
                                    elseif name == "yorigin"
                                        origin_y = parse(Float64, content(fprop))
                                    elseif name == "xscale"
                                        scale = parse(Float64, content(fprop))
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    mass_positions = Dict{String, Vector{Tuple{Float64, Float64}}}()
    
    for prop in child_elements(xroot)
        if attribute(prop, "name") == "tracks"
            for item in child_elements(prop)
                for obj in child_elements(item)
                    classname = attribute(obj, "class")
                    if classname !== nothing && contains(classname, "PointMass")
                        mass_name = ""
                        positions = Tuple{Float64, Float64}[]
                        
                        for mprop in child_elements(obj)
                            if attribute(mprop, "name") == "name"
                                mass_name = content(mprop)
                            elseif attribute(mprop, "name") == "framedata"
                                for frame in child_elements(mprop)
                                    x, y = 0.0, 0.0
                                    for fobj in child_elements(frame)
                                        for fprop in child_elements(fobj)
                                            if attribute(fprop, "name") == "x"
                                                x = parse(Float64, content(fprop))
                                            elseif attribute(fprop, "name") == "y"
                                                y = parse(Float64, content(fprop))
                                            end
                                        end
                                    end
                                    x_real = (x - origin_x) / scale
                                    y_real = (y - origin_y) / scale
                                    push!(positions, (x_real, y_real))
                                end
                            end
                        end
                        
                        if !isempty(mass_name) && !isempty(positions)
                            mass_positions[mass_name] = positions
                        end
                    end
                end
            end
        end
    end
    
    free(xdoc)
    return mass_positions, origin_x, origin_y, scale
end

function unwrap_angles(angles)
    unwrapped = mod.(angles, 2π)
    
    for i in 2:length(unwrapped)
        diff = unwrapped[i] - unwrapped[i-1]
        if diff > π
            unwrapped[i:end] .-= 2π
        elseif diff < -π
            unwrapped[i:end] .+= 2π
        end
    end
    return unwrapped
end

function calculate_angles_from_positions(mass_A_positions, mass_B_positions)
    n = min(length(mass_A_positions), length(mass_B_positions))
    θ1_raw = zeros(n)
    θ2_raw = zeros(n)
    
    for i in 1:n
        x1, y1 = mass_A_positions[i]
        x2, y2 = mass_B_positions[i]
        
        θ1_raw[i] = atan(x1, y1)
        
        dx = x2 - x1
        dy = y2 - y1
        θ2_raw[i] = atan(dx, dy)
    end
    
    θ1 = unwrap_angles(θ1_raw)
    θ2 = unwrap_angles(θ2_raw)
    
    return θ1, θ2
end

function calculate_initial_velocities(θ1_data, θ2_data, dt)
    ω1_0 = (θ1_data[3] - θ1_data[1]) / (2 * dt)
    ω2_0 = (θ2_data[3] - θ2_data[1]) / (2 * dt)
    return ω1_0, ω2_0
end

# =============================================================================
# Fonctions physiques
# =============================================================================

function position(θ1, L1, θ2, L2)
    r1 = [L1*sin(θ1), L1*cos(θ1)]
    r2 = [r1[1] + L2*sin(θ2), r1[2] + L2*cos(θ2)]
    return r1, r2
end

function equations_double_pendulum(du, u, p, t)
    θ1, ω1, θ2, ω2 = u
    m1, m2, L1, L2, g = p
    
    δ = θ2 - θ1
    cos_δ = cos(δ)
    sin_δ = sin(δ)
    
    den = m1 + m2 * sin_δ^2
    
    du[1] = ω1
    du[2] = (m2 * sin_δ * (L1 * ω1^2 * cos_δ + L2 * ω2^2) +
             m2 * g * sin(θ2) * cos_δ -
             (m1 + m2) * g * sin(θ1)) / (L1 * den)
    
    du[3] = ω2
    du[4] = (-m2 * L2 * ω2^2 * sin_δ * cos_δ -
             (m1 + m2) * L1 * ω1^2 * sin_δ +
             (m1 + m2) * g * sin(θ1) * cos_δ -
             (m1 + m2) * g * sin(θ2)) / (L2 * den)
end

function simulate_pendulum(m1, m2, L1, L2, g, θ1_0, ω1_0, θ2_0, ω2_0, tspan, saveat)
    u0 = [θ1_0, ω1_0, θ2_0, ω2_0]
    p = [m1, m2, L1, L2, g]
    
    prob = ODEProblem(equations_double_pendulum, u0, tspan, p)
    sol = solve(prob, Tsit5(), saveat=saveat)
    
    return sol
end

# =============================================================================
# Paramètres physiques
# =============================================================================

L1 = 91.74 / 1000 # [m]
L2 = 69.33 / 1000 # [m]
g = 9.81 # [m/s²]

# =============================================================================
# Charger les données Tracker
# =============================================================================

println("Chargement des données Tracker...")
trk_path = joinpath(@__DIR__, "assets", "First_Video_2s.trk")
mass_positions, origin_x, origin_y, scale = load_tracker_positions(trk_path)

dt = 0.01
n_frames = 200
t_data = collect(0:dt:(n_frames-1)*dt)

if haskey(mass_positions, "mass A") && haskey(mass_positions, "mass B")
    mass_A = mass_positions["mass A"]
    mass_B = mass_positions["mass B"]
    
    n_actual = min(length(mass_A), length(mass_B), n_frames)
    t_data = t_data[1:n_actual]
    
    θ1_data, θ2_data = calculate_angles_from_positions(mass_A[1:n_actual], mass_B[1:n_actual])
    
    θ1_0 = θ1_data[1]
    θ2_0 = θ2_data[1]
    
    ω1_init, ω2_init = calculate_initial_velocities(θ1_data, θ2_data, dt)
    
    println("Conditions initiales estimées:")
    println("  θ1_0 = $(round(rad2deg(θ1_0), digits=1))°, θ2_0 = $(round(rad2deg(θ2_0), digits=1))°")
    println("  ω1_0 = $(round(rad2deg(ω1_init), digits=1))°/s, ω2_0 = $(round(rad2deg(ω2_init), digits=1))°/s")
else
    error("Masses 'mass A' et 'mass B' non trouvées!")
end

# =============================================================================
# Optimisation
# =============================================================================

t_optim_end = 0.5
idx_optim = findlast(t -> t <= t_optim_end, t_data)
t_optim = t_data[1:idx_optim]
θ1_optim = θ1_data[1:idx_optim]
θ2_optim = θ2_data[1:idx_optim]

function objective_full(params)
    m1, m2, ω1_0, ω2_0 = params
    
    sol = simulate_pendulum(m1, m2, L1, L2, g, θ1_0, ω1_0, θ2_0, ω2_0, 
                           (t_optim[1], t_optim[end]), t_optim)
    
    error = 0.0
    for i in 1:length(t_optim)
        error += (sol[1, i] - θ1_optim[i])^2 + (sol[3, i] - θ2_optim[i])^2
    end
    return error
end

println("\nOptimisation de [m1, m2, ω1_0, ω2_0] sur $(t_optim_end)s...")

m1_init = 0.004
m2_init = 0.003

lower = [0.001, 0.001, ω1_init - 3.0, ω2_init - 3.0]
upper = [0.01,   0.01,   ω1_init + 3.0, ω2_init + 3.0]
initial = [m1_init, m2_init, ω1_init, ω2_init]

result = optimize(
    objective_full,
    lower,
    upper,
    initial,
    Fminbox(NelderMead()),
    Optim.Options(iterations=10000, show_trace=false)
)

m1_opt = Optim.minimizer(result)[1]
m2_opt = Optim.minimizer(result)[2]
ω1_opt = Optim.minimizer(result)[3]
ω2_opt = Optim.minimizer(result)[4]

println("\nRésultats de l'optimisation:")
println("  m1 = $(round(m1_opt*1000, digits=1)) g, m2 = $(round(m2_opt*1000, digits=1)) g (ratio (m1+m2)/m1 = $(round((m1_opt+m2_opt)/m1_opt, digits=2)))")
println("  ω1_0 = $(round(rad2deg(ω1_opt), digits=1))°/s, ω2_0 = $(round(rad2deg(ω2_opt), digits=1))°/s")
println("  Erreur finale = $(round(Optim.minimum(result), digits=6))")

# =============================================================================
# Simulation finale
# =============================================================================

t_extension = 2.0
t_final = t_data[end] + t_extension
t_sim_extended = collect(0:dt:t_final)

println("\nSimulation sur $(round(t_final, digits=1))s (données réelles: $(round(t_data[end], digits=1))s + extrapolation: $(t_extension)s)")

sol_opt = simulate_pendulum(m1_opt, m2_opt, L1, L2, g, θ1_0, ω1_opt, θ2_0, ω2_opt,
                            (0.0, t_final), t_sim_extended)

global θ1_sim_full = sol_opt[1, :]
global θ2_sim_full = sol_opt[3, :]
global t_sim_full = t_sim_extended

n_data = length(θ1_data)

for i in 1:n_data
    diff = θ2_sim_full[i] - θ2_data[i]
    n_tours = round(diff / (2π))
    θ2_sim_full[i] -= n_tours * 2π
end

global θ1_sim = θ1_sim_full[1:n_data]
global θ2_sim = θ2_sim_full[1:n_data]

sol_opt_stats = simulate_pendulum(m1_opt, m2_opt, L1, L2, g, θ1_0, ω1_opt, θ2_0, ω2_opt,
                                   (t_data[1], t_data[end]), t_data)

m_opt = [m1_opt, m2_opt]
u0 = [θ1_0, ω1_opt, θ2_0, ω2_opt]
p_opt = [m1_opt, m2_opt, L1, L2, g]
prob_opt = ODEProblem(equations_double_pendulum, u0, (t_data[1], t_data[end]), p_opt)
sol_opt = sol_opt_stats

println("\n✓ Prêt pour statistiques.jl")