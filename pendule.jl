include("imports.jl")

# =============================================================================
# Chargement des données Tracker (fichier .trk)
# =============================================================================

function load_tracker_positions(filepath::String)
    """Charge les positions des masses depuis le fichier .trk"""
    xdoc = parse_file(filepath)
    xroot = root(xdoc)
    
    # Paramètres de calibration
    origin_x = 0.0
    origin_y = 0.0
    scale = 1.0
    
    # Parcourir pour trouver ImageCoordSystem
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
    
    # Extraire les données des PointMass
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
                                    # Convertir en coordonnées réelles (mètres)
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
    """
    Corrige les sauts de 2π dans une série d'angles pour obtenir une trajectoire continue.
    """
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
        
        # θ1 : angle du premier bras depuis la verticale
        θ1_raw[i] = atan(x1, y1)
        
        # θ2 : angle du DEUXIÈME bras PAR RAPPORT AU PREMIER (angle relatif)
        dx = x2 - x1
        dy = y2 - y1
        θ2_raw[i] = atan(dx, dy)
    end
    
    θ1 = unwrap_angles(θ1_raw)
    θ2 = unwrap_angles(θ2_raw)
    
    return θ1, θ2
end

function calculate_initial_velocities(θ1_data, θ2_data, dt)
    """
    Calcule les vitesses angulaires initiales par différences finies centrées
    """
    # Utiliser les 3 premiers points pour une meilleure estimation
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

function equations_double_pendulum!(du, u, p, t)
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
    
    prob = ODEProblem(equations_double_pendulum!, u0, tspan, p)
    sol = solve(prob, Tsit5(), saveat=saveat)
    
    return sol
end

# =============================================================================
# Paramètres physiques
# =============================================================================

L1 = 91.74 / 1000  # m
L2 = 69.33 / 1000  # m
g = 9.81           # m/s²

# =============================================================================
# Charger les données Tracker
# =============================================================================

println("Chargement des données Tracker...")
trk_path = joinpath(@__DIR__, "assets", "First_Video_2s.trk")
mass_positions, origin_x, origin_y, scale = load_tracker_positions(trk_path)

println("  Origine: ($origin_x, $origin_y) pixels")
println("  Échelle: $scale pixels/m")
println("  Masses trouvées: $(keys(mass_positions))")

dt = 0.01
n_frames = 200
t_data = collect(0:dt:(n_frames-1)*dt)

if haskey(mass_positions, "mass A") && haskey(mass_positions, "mass B")
    mass_A = mass_positions["mass A"]
    mass_B = mass_positions["mass B"]
    
    n_actual = min(length(mass_A), length(mass_B), n_frames)
    t_data = t_data[1:n_actual]
    
    θ1_data, θ2_data = calculate_angles_from_positions(mass_A[1:n_actual], mass_B[1:n_actual])
    # Après le chargement des données, avant l'optimisation
    println("\n=== DIAGNOSTIC GÉOMÉTRIQUE ===")
    for frame in [1, 10, 50]
        x1, y1 = mass_A[frame]
        x2, y2 = mass_B[frame]
        
        # Calculer les deux versions
        θ1_abs = atan(x1, y1)
        θ2_abs = atan(x2, y2)          # Angle absolu
        θ2_rel = atan(x2 - x1, y2 - y1) # Angle relatif
        
        println("\nFrame $frame:")
        println("  Position A: ($(round(x1,digits=4)), $(round(y1,digits=4)))")
        println("  Position B: ($(round(x2,digits=4)), $(round(y2,digits=4)))")
        println("  θ1 = $(round(rad2deg(θ1_abs), digits=1))°")
        println("  θ2 absolu = $(round(rad2deg(θ2_abs), digits=1))°")
        println("  θ2 relatif = $(round(rad2deg(θ2_rel), digits=1))°")
        println("  Différence θ2_abs - θ1 = $(round(rad2deg(θ2_abs - θ1_abs), digits=1))°")
        println("  → θ2_relatif devrait ≈ (θ2_abs - θ1)")
    end
    
    θ1_0 = θ1_data[1]
    θ2_0 = θ2_data[1]
    
    # Calculer vitesses initiales estimées
    ω1_init, ω2_init = calculate_initial_velocities(θ1_data, θ2_data, dt)
    
    println("\nConditions initiales estimées:")
    println("  θ1_0 = $(rad2deg(θ1_0))°")
    println("  θ2_0 = $(rad2deg(θ2_0))°")
    println("  ω1_0 = $(rad2deg(ω1_init))°/s")
    println("  ω2_0 = $(rad2deg(ω2_init))°/s")
    println("  Nombre de frames: $n_actual")
else
    error("Masses 'mass A' et 'mass B' non trouvées dans le fichier Tracker!")
end

# =============================================================================
# Optimisation
# =============================================================================

t_optim_end = 0.5  # Optimiser sur 0.5s (fenêtre raisonnable avant divergence chaotique)
idx_optim = findlast(t -> t <= t_optim_end, t_data)
t_optim = t_data[1:idx_optim]
θ1_optim = θ1_data[1:idx_optim]
θ2_optim = θ2_data[1:idx_optim]

function objective_full(params)
    """
    Optimise m1, m2, ω1_0, ω2_0
    params = [m1, m2, ω1_0, ω2_0]
    """
    m1 = params[1]
    m2 = params[2]
    ω1_0 = params[3]
    ω2_0 = params[4]
    
    sol = simulate_pendulum(m1, m2, L1, L2, g, θ1_0, ω1_0, θ2_0, ω2_0, 
                           (t_optim[1], t_optim[end]), t_optim)
    
    error = 0.0
    for i in 1:length(t_optim)
        # Poids égaux pour θ1 et θ2
        error += (sol[1, i] - θ1_optim[i])^2 + (sol[3, i] - θ2_optim[i])^2
    end
    return error
end

println("\n" * "="^60)
println("OPTIMISATION DES PARAMÈTRES")
println("="^60)
println("\nOptimisation de [m1, m2, ω1_0, ω2_0] sur $(t_optim_end)s...")

# Valeurs initiales
m1_init = 0.05  # kg
m2_init = 0.03  # kg

# Bornes de recherche
lower = [0.001, 0.001, ω1_init - 3.0, ω2_init - 3.0]
upper = [0.2,   0.2,   ω1_init + 3.0, ω2_init + 3.0]
initial = [m1_init, m2_init, ω1_init, ω2_init]

result = optimize(
    objective_full,
    lower,
    upper,
    initial,
    Fminbox(NelderMead()),
    Optim.Options(iterations=10000, show_trace=false)
)

# Résultats optimaux
m1_opt = Optim.minimizer(result)[1]
m2_opt = Optim.minimizer(result)[2]
ω1_opt = Optim.minimizer(result)[3]
ω2_opt = Optim.minimizer(result)[4]

println("\n" * "="^60)
println("RÉSULTATS DE L'OPTIMISATION")
println("="^60)
println("\nMasses optimales:")
println("  m1 = $(round(m1_opt*1000, digits=2)) g")
println("  m2 = $(round(m2_opt*1000, digits=2)) g")
println("  Ratio m2/m1 = $(round(m2_opt/m1_opt, digits=3))")

println("\nVitesses angulaires initiales optimales:")
println("  ω1_0 = $(round(rad2deg(ω1_opt), digits=2))°/s (estimé: $(round(rad2deg(ω1_init), digits=2))°/s)")
println("  ω2_0 = $(round(rad2deg(ω2_opt), digits=2))°/s (estimé: $(round(rad2deg(ω2_init), digits=2))°/s)")

println("\nDifférence avec estimation initiale:")
println("  Δω1 = $(round(rad2deg(ω1_opt - ω1_init), digits=2))°/s")
println("  Δω2 = $(round(rad2deg(ω2_opt - ω2_init), digits=2))°/s")

println("\nErreur finale = $(round(Optim.minimum(result), digits=6))")
println("Convergence: $(Optim.converged(result))")

# =============================================================================
# Simulation finale
# =============================================================================

println("\n" * "="^60)
println("SIMULATION FINALE")
println("="^60)

# Simulation avec paramètres optimisés
sol_opt = simulate_pendulum(m1_opt, m2_opt, L1, L2, g, θ1_0, ω1_opt, θ2_0, ω2_opt,
                            (t_data[1], t_data[end]), t_data)

θ1_sim = sol_opt[1, :]
θ2_sim = sol_opt[3, :]
# Corriger les tours multiples pour θ2
println("\n=== CORRECTION DES TOURS MULTIPLES ===")
println("Avant correction:")
println("  θ2_sim[end] = $(rad2deg(θ2_sim[end]))°")
println("  θ2_data[end] = $(rad2deg(θ2_data[end]))°")

# Ajuster pour matcher le nombre de tours
for i in 1:length(θ2_sim)
    diff = θ2_sim[i] - θ2_data[i]
    n_tours = round(diff / (2π))
    θ2_sim[i] -= n_tours * 2π
end

println("\nAprès correction:")
println("  θ2_sim[end] = $(rad2deg(θ2_sim[end]))°")
println("  θ2_data[end] = $(rad2deg(θ2_data[end]))°")
println("  Erreur finale = $(rad2deg(abs(θ2_sim[end] - θ2_data[end])))°")

# Variables pour compatibilité avec statistiques.jl
m_opt = [m1_opt, m2_opt]
u0 = [θ1_0, ω1_opt, θ2_0, ω2_opt]
p_opt = [m1_opt, m2_opt, L1, L2, g]
prob_opt = ODEProblem(equations_double_pendulum!, u0, (t_data[1], t_data[end]), p_opt)

println("\nOptimisation effectuée sur $(round(t_optim_end, digits=2))s ($idx_optim frames)")
println("Simulation finale sur $(round(t_data[end], digits=2))s ($(length(t_data)) frames)")
println("\nExécutez statistiques.jl pour les analyses détaillées.")
println("="^60)