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
    Convertit d'abord en base 0-2π (360°), puis déroule les discontinuités.
    """
    # D'abord convertir en base 0-2π
    unwrapped = mod.(angles, 2π)
    
    # Puis dérouler les sauts
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
    """
    Calcule les angles θ1 et θ2 à partir des positions des masses.
    
    Convention Tracker (confirmée):
    - x positif = vers la droite  
    - y positif = vers le HAUT
    - θ = 0 quand le pendule est en bas (y négatif)
    
    Convention physique standard du pendule double:
    - θ = 0 quand le pendule pointe vers le bas
    - x = L*sin(θ), y = -L*cos(θ)
    
    Dans Tracker: pendule en bas → y négatif, pendule en haut → y positif
    Pour atan2(x, y): on obtient θ=0 quand x=0 et y>0
    
    Donc pour avoir θ=0 quand le pendule est en bas (y<0 dans Tracker):
    θ = atan(x, y) donne directement le bon angle car:
    - Pendule en bas: y < 0 → atan(0, y<0) = π (ou -π), NON!
    
    En fait, si y positif = haut et on veut θ=0 en bas:
    θ = atan(x, y) avec y négatif en bas → atan(0, -L) = π
    
    Mais Tracker dit θ=0 en bas... donc Tracker utilise atan(x, -y)!
    Non attendez - revoyons: si Tracker montre θ=0 en bas, c'est SA convention.
    
    Pour notre modèle physique où y négatif = bas:
    y_physique = y_tracker (même convention!)
    θ = atan(x, y) directement
    
    Note: Les angles sont "unwrapped" pour éviter les sauts à ±180°
    """
    n = min(length(mass_A_positions), length(mass_B_positions))
    θ1_raw = zeros(n)
    θ2_raw = zeros(n)
    
    for i in 1:n
        x1, y1 = mass_A_positions[i]
        x2, y2 = mass_B_positions[i]
        
        # Convention Tracker: y positif = haut, θ=0 en haut
        # atan(x, y) donne l'angle depuis l'axe +y (antihoraire = positif)
        θ1_raw[i] = atan(x1, y1)
        
        # θ2: angle du deuxième pendule par rapport à la verticale
        dx = x2 - x1
        dy = y2 - y1
        θ2_raw[i] = atan(dx, dy)
    end
    
    # Unwrap pour avoir des angles continus (corrige les sauts à ±180°)
    θ1 = unwrap_angles(θ1_raw)
    θ2 = unwrap_angles(θ2_raw)
    
    return θ1, θ2
end

# Charger les données réelles
println("Chargement des données Tracker...")
trk_path = joinpath(@__DIR__, "assets", "First_Video_2s.trk")
mass_positions, origin_x, origin_y, scale = load_tracker_positions(trk_path)

println("  Origine: ($origin_x, $origin_y) pixels")
println("  Échelle: $scale pixels/m")
println("  Masses trouvées: $(keys(mass_positions))")

# Paramètres temporels (100 fps, 2 secondes = 200 frames)
dt = 0.01  # 10ms entre chaque frame
n_frames = 200
t_data = collect(0:dt:(n_frames-1)*dt)

# =============================================================================
# Fonctions utilitaires
# =============================================================================

function degrees_to_radians(degrees)
    return degrees * (π / 180)
end

function Lagrangien(θ1_0, θ2_0, m1, m2, g, r1, r2, L1, L2)
    lagrange = m1*(L1*θ1_0)^2 /2 + m2*((L1*θ1_0)^2 + (L2*θ2_0)^2 + 2*L1*L2*θ1_0*θ2_0*cos(θ1_0 - θ2_0)) + m1*g*L1*cos(θ1_0) + m2*g*(L1*cos(θ1_0) + L2*cos(θ2_0))/2
    return lagrange
end

function position(θ1_0, L1, θ2_0, L2)
    r1 = [L1*sin(θ1_0), L1*cos(θ1_0)]
    r2 = [r1[1] + L2*sin(θ2_0), r1[2] + L2*cos(θ2_0)]
    return r1, r2
end

function equations_double_pendulum!(du, u, p, t)
    θ1, ω1, θ2, ω2 = u
    m1, m2, L1, L2, g = p
    
    δ = θ2 - θ1
    cos_δ = cos(δ)
    sin_δ = sin(δ)
    
    # Dénominateur commun (ne s'annule jamais car m1 + m2*sin²(δ) > 0)
    den = m1 + m2 * sin_δ^2
    
    # Accélérations angulaires (équations de Lagrange correctes)
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

function simulate_pendulum(masses, L1, L2, g, θ1_0, θ2_0, tspan, saveat)
    m1, m2 = masses
    u0 = [θ1_0, 0.0, θ2_0, 0.0]  # conditions initiales (angles, vitesses angulaires)
    p = [m1, m2, L1, L2, g]  # L1 et L2 déjà en mètres
    
    prob = ODEProblem(equations_double_pendulum!, u0, tspan, p)
    sol = solve(prob, Tsit5(), saveat=saveat)
    
    return sol
end

function cost_function(masses, L1, L2, g, θ1_0, θ2_0, t_data, θ1_data, θ2_data)
    sol = simulate_pendulum(masses, L1, L2, g, θ1_0, θ2_0, (t_data[1], t_data[end]), t_data)
    
    # Erreur quadratique (somme des carrés des résidus)
    error = 0.0
    for i in 1:length(t_data)
        error += (sol[1, i] - θ1_data[i])^2 + (sol[3, i] - θ2_data[i])^2
    end
    
    return error
end

# =============================================================================
# Paramètres physiques
# =============================================================================

# ATTENTION: Les longueurs sont en mm, on les convertit en m
L1 = 91.74 / 1000  # Longueur du premier pendule (m) - était 91.74mm
L2 = 69.33 / 1000  # Longueur du second pendule (m) - était 69.33mm

g = 9.81    # Accélération due à la gravité (m/s²)

# =============================================================================
# Charger les vraies données depuis Tracker
# =============================================================================

# Extraire les positions et calculer les angles
if haskey(mass_positions, "mass A") && haskey(mass_positions, "mass B")
    mass_A = mass_positions["mass A"]
    mass_B = mass_positions["mass B"]
    
    # Ajuster le nombre de frames si nécessaire
    n_actual = min(length(mass_A), length(mass_B), n_frames)
    t_data = t_data[1:n_actual]
    
    # Calculer les angles depuis les positions
    θ1_data, θ2_data = calculate_angles_from_positions(mass_A[1:n_actual], mass_B[1:n_actual])
    
    # Angles initiaux (directement depuis les données)
    θ1_0 = θ1_data[1]
    θ2_0 = θ2_data[1]
    
    println("\nAngles initiaux:")
    println("  θ1_0 = $(rad2deg(θ1_0))°")
    println("  θ2_0 = $(rad2deg(θ2_0))°")
    println("  Nombre de frames: $n_actual")
else
    error("Masses 'mass A' et 'mass B' non trouvées dans le fichier Tracker!")
end

# Utiliser Optim.jl directement pour NelderMead (plus simple et sans problème de gradients)
# Optimiser sur une fenêtre plus large pour mieux capturer la dynamique
t_optim_end = 1.0  # Augmenté de 0.02s à 1.0s pour mieux coller aux données
idx_optim = findlast(t -> t <= t_optim_end, t_data)
t_optim = t_data[1:idx_optim]
θ1_optim = θ1_data[1:idx_optim]
θ2_optim = θ2_data[1:idx_optim]

function objective(masses)
    # Pondération beaucoup plus forte sur θ2 pour réduire son erreur
    sol = simulate_pendulum(masses, L1, L2, g, θ1_0, θ2_0, (t_optim[1], t_optim[end]), t_optim)
    error = 0.0
    for i in 1:length(t_optim)
        error += (sol[1, i] - θ1_optim[i])^2 + 5.0 * (sol[3, i] - θ2_optim[i])^2
    end
    return error
end

m_init = [0.05, 0.03]  # kg - valeurs initiales raisonnables

println("\nOptimisation des masses...")
# Utiliser optimize de Optim.jl directement
result = optimize(objective, [0.001, 0.001], [1.0, 1.0], m_init, Fminbox(NelderMead()))

m_opt = Optim.minimizer(result)

println("\nMasses optimales:")
println("  m1 = $(m_opt[1]*1000) g")
println("  m2 = $(m_opt[2]*1000) g")
println("  Ratio m2/m1 = $(m_opt[2]/m_opt[1])")
println("  Erreur finale = $(Optim.minimum(result))")

# Simulation finale avec les masses optimisées
u0 = [θ1_0, 0.0, θ2_0, 0.0]
p_opt = [m_opt[1], m_opt[2], L1, L2, g]
prob_opt = ODEProblem(equations_double_pendulum!, u0, (t_data[1], t_data[end]), p_opt)
sol_opt = solve(prob_opt, Tsit5(), saveat=t_data)

θ1_sim = sol_opt[1, :]
θ2_sim = sol_opt[3, :]

println("\nOptimisation sur les $(t_optim_end*1000) premières ms ($idx_optim frames)")
println("Simulation terminée sur $(t_data[end])s ($(length(t_data)) frames).")
println("Exécutez statistiques.jl pour les analyses détaillées.")