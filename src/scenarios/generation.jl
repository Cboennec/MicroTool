#=
    Scenario generation methods
=#
abstract type AbstractScenariosGenerator end

mutable struct MarkovChain
    states
    transition_matrices
end

# Compute markov chain
#########################################################################################################################################################################################
# On construit ici les états clusterisé ainsi que l'affectation des data points originaux aux clusters (voir  compute_markov_states(data, nstate::Int64, algo::String) )
# Puis en utilisant l'affectation des data points des series originales aux différent clusters (stocké dans sequence) on construit les probabilité de transition (voir compute_transition_matrices_between_hours(sequence, nstate::Int64) )
# Enfin on a tout les éléments pour renvoyer une instance de la structure MarkovChain
#########################################################################################################################################################################################
function MarkovChain(data, nstate::Int64, algo::String)
    # Parameters
    nk = size(data,1) # number of data series (3)
    nh = 24 # number of hours (24)
    nm = size(data[1],1) # number of month (12)

    # Pre-allocate
    states = [zeros(nh, nk, nstate) for m in 1:nm]
    transition_matrices = [zeros(nstate, nstate, nh - 1) for m in 1:nm]


    # For each month, we compute the markov states and transition matrices
    for m in 1:nm
        # For each hour, we extract markov states using k-means from aggregated data to keep the synchronicity
        states[m], sequences = compute_markov_states([data[k][m] for k in 1:nk], nstate, algo)
        # We compute the transition matrices between consecutive hours using the sequences
        transition_matrices[m] = compute_transition_matrices_between_hours(sequences, nstate)
    end

    return MarkovChain(states, transition_matrices)
end
# Compute states using clustering algorithm
#########################################################################################################################################################################################
# Ici la donnée d'entrée data contient 1 matrice par serie de données. Les matrices sont de tille 24*N (24h fois N, le nombre de jour à disposition pour ce mois ce type de jour (wk/wkd) sur X années )
# D'abord on normalise les données
# Puis pour chaque heure on va créer N data points chacun un état (pv, ld_E, ld_H)
# On affecte chaque data point à un cluster parmis nstate clusters
# On dénormalise et on stock dans states pour chaque heure, les valeurs (pv, ld_E, ld_H) de chaque nstate états.
# Enfin on stock dans sequence la suite d'affectation des N data points aux clusters
#
# Note : kmeans et kmedoids viennent du package Clustering
# Note 2 : Ici la matrice de séquence a 25 lignes pour pouvoir faire que la première et la dernière soient égales et ainsi former une forme de boucle.
# C'est un parti pris qui pourra être modifié
#########################################################################################################################################################################################
function compute_markov_states(data, nstate::Int64, algo::String)
    # Parameters
    nk = size(data,1)  # number of data series (3)
    nh = size(data[1], 1) # number of hours (24)
    nd = size(data[1], 2) # number of available days for this month and type of day called N ou N_obs

    # Pre-allocate
    states = zeros(nh, nk, nstate)
    sequences = zeros(Int64, nh+1, nd)

    # Normalization
    data_n = replace!.(data ./ maximum.(data), NaN => 0.)

    # Clustering
    for h in 1:nh
        if algo == "kmeans"
            # Aggregate normalized data
            data_agg = permutedims(hcat([data_n[k][h,:] for k in 1:nk]...))
            # Clustering
            clusters = kmeans(data_agg, nstate)
            # Store the denormalized states and sequences
            states[h,:,:] = clusters.centers .* maximum.(data)
            sequences[h,:] = permutedims(clusters.assignments)
        elseif algo == "kmedoids"
            # Aggregate normalized data
            data_agg = permutedims(hcat([data_n[k][h,:] for k in 1:nk]...))
            # Compute euclidean distance matrix
            dist = pairwise(Euclidean(), data_agg, dims = 2 )
            # Clustering
            clusters = kmedoids(dist, nstate)
            # Store the states and sequences
            states[h,:,:] = data_agg[:, clusters.medoids] .* maximum.(data)
            sequences[h,:] = permutedims(clusters.assignments)
        end
    end
    # The sequence of the last hour is equal to the first one
    sequences[end, :] = sequences[1, :]

    return states, sequences
end
# Compute transition matrix between consecutive hours
#########################################################################################################################################################################################
# Ici on cherche a calculer la probabilité de passer de chaque nstate cluster de l'heure h à chaque nstate cluster de l'heure h+1
# Pour chaque transition (pour chaque heure)
#   On incrémente chaque transition (i,j,h) du cluster i à l'heure h vers le cluster j à l'heure h+1
#Puis on divise les nombre obtenus par le nombre de data point qui appartenait à ce cluster ( = somme des points partant d'un cluster = somme sur la ligne du clisuter concerné)
#Pour une explication plus profonde et illustrée voir pdf (probablement figure 7: Defining transition matrix probabilities)
#########################################################################################################################################################################################
function compute_transition_matrices_between_hours(sequence, nstate::Int64)
    # Parameters
    nh = size(sequence, 1) # number of step in the chain

    # Pre-allocate
    transition_matrix = zeros(nstate, nstate, nh - 1)

    for h in 1:nh-1
        for (i,j) in zip(sequence[h, :], sequence[h+1, :])
            transition_matrix[i, j, h] += 1
        end

        # Convert into probabilities
        for i in 1:size(transition_matrix,1)

            sum_row = sum(transition_matrix[i, :, h])

            if sum_row > 0 # to avoid zero division...
                transition_matrix[i, :, h] = transition_matrix[i, :, h] ./ sum_row
            elseif sum_row == 0
                transition_matrix[i, i, h] = 1 # to avoid bug in Categorical...
            else
                print("Error! The sum couldn't be negative")
                return
            end
        end
    end

    return transition_matrix
end

# Clustering data by week for each month
#########################################################################################################################################################################################
#Separe les données en 12 vecteur (1 par mois) contenant les données de puissance des jours de semaine ( plus exactement qui ne sont pas du weekend)
#On restructure ces données en matrice 24*N afin de les regrouper par jour où N est le nombre de jour de semaine pour ce mois (sur X années de données, ici 3 années)

#Note : isweekend est dans scenarios/utils.jl et renvoi True pour un (samedi || dimanche) et False sinon
#########################################################################################################################################################################################
function clustering_month_week(data, t; flag::String="")
    # Parameters
    nm = 12 # number of month in one year

    # Pre_allocate
    data_cluster= Vector{Array{Float64,2}}(undef, nm)

    for m in 1:nm
        # We cluster week days using the "not" operator
        if flag == "odd" # odd indices
            data_cluster[m] = reshape(data[(Dates.month.(t) .== m) .& .~isweekend(t)], 24, :)[:, isodd.(1:end)]
        elseif flag == "even" # even indices
            data_cluster[m] = reshape(data[(Dates.month.(t) .== m) .& .~isweekend(t)], 24, :)[:, iseven.(1:end)]
        else
            data_cluster[m] = reshape(data[(Dates.month.(t) .== m) .& .~isweekend(t)], 24, :)
        end
    end

    return data_cluster
end
# Clustering data by weekend for each month
#########################################################################################################################################################################################
#Separe les données en 12 vecteur (1 par mois) contenant les données de puissance des jours de weekend
#On restructure ces données en matrice 24*N afin de les regrouper par jour où N est le nombre de jour de weekend pour ce mois (sur X années de données, ici 3 années)

#Note : isweekend est dans scenarios/utils.jl et renvoi True pour un (samedi || dimanche) et False sinon
#########################################################################################################################################################################################

function clustering_month_weekend(data, t; flag::String="")
    # Parameters
    nm = 12 # number of month in one year

    # Pre_allocate
    data_cluster = Vector{Array{Float64,2}}(undef, nm)

    for m in 1:nm
        # We cluster weekend days
        if flag == "odd" # odd indices
            data_cluster[m] = reshape(data[(Dates.month.(t) .== m) .& isweekend(t)], 24, :)[:, isodd.(1:end)]
        elseif flag == "even" # even indices
            data_cluster[m] = reshape(data[(Dates.month.(t) .== m) .& isweekend(t)], 24, :)[:, iseven.(1:end)]
        else
            data_cluster[m] = reshape(data[(Dates.month.(t) .== m) .& isweekend(t)], 24, :)
        end
    end

    return data_cluster
end
# Clustering data by month
function clustering_month(data, t; flag::String="")
    # Parameters
    nm = 12 # number of month in one year

    # Pre_allocate
    data_cluster = Vector{Array{Float64,2}}(undef, nm)

    for m in 1:nm
        # We cluster weekend days
        if flag == "odd" # odd indices
            data_cluster[m] = reshape(data[(Dates.month.(t) .== m)], 24, :)[:, isodd.(1:end)]
        elseif flag == "even" # even indices
            data_cluster[m] = reshape(data[(Dates.month.(t) .== m)], 24, :)[:, iseven.(1:end)]
        else
            data_cluster[m] = reshape(data[(Dates.month.(t) .== m)], 24, :)
        end
    end

    return data_cluster
end

mutable struct MarkovGenerator <: AbstractScenariosGenerator
    nstate::Int64 # Nombre d'états et par conséquent nombre de cluster par heure par mois.
    algo::String # Algorithme utilisé pour le clustering des états.
    markovchains::NamedTuple{(:wk, :wkd), Tuple{MarkovChain, MarkovChain}}

    MarkovGenerator(; nstate = 10, algo = "kmeans") = new(nstate, algo)
end

########################################################################################################################################################################################
# On recupère dans wk et wkd pour chaque série de données une matrice par mois 24*N, N étant le nombre de jours à disposition pour ce mois en semaine/weekend. On donc pour 24h dans chaue mois N data points
# nstate désigne le nombre de cluster qui vont être constitués, on vérifie donc que sa valeure est pertinente  (on ne peut pas faire plus de cluster qu'il n'y a de data point)
# Ici generator à un attribut markovchains qui lui même est un tuple contenant 2 structures markovchain labelisées (:wk et :wkd)
# Voir la fonction markovchain qui construit et renvoi une structure de type markovchain avec en entrée les série de donnée par heure par mois par jour de semaine/weekend, un nombre de cluster par étape, un algorithme de clusterisation.
########################################################################################################################################################################################

function initialize_generator!(generator::MarkovGenerator, data...)
    # TODO: generic function without .power and .t or specify the input type

    # Clustering by week and weekend for each month for each data series (pv, ld_E, ld_H)
    wk = [clustering_month_week(d.power, d.t) for d in data]
    wkd = [clustering_month_weekend(d.power, d.t) for d in data]

    # The number of states must be greater than 2 and lower than the minimum number of week and weekend days
    min_wk, min_wkd =  minimum([size(wk[k][m])[2] for k in 1:size(wk,1), m in 1:12]), minimum([size(wkd[k][m])[2] for k in 1:size(wk,1), m in 1:12])
    1 < generator.nstate < min(min_wk, min_wkd) ? nothing : generator.nstate = min(min_wk, min_wkd) - 1

    # Compute markov chains

    generator.markovchains = (wk = MarkovChain(wk, generator.nstate, generator.algo),
                             wkd = MarkovChain(wkd, generator.nstate, generator.algo))

    return generator
end

# Generate scenario from markov chains
######################################################################################################################################################################################
# Param d'entrée :
#  -generator : voir la structure de donnée associée et la fonction initialize_generator!
#  -s0 : "état initial"
#  -t0 : la date de départ
#  -nstep : longueur de la séquence à générer par an (ici une année en heure par heure = 8760 étapes)
#  -ny : nombre d'année des scénario à générer
#  -ns : nombre de scénario a générer

# On commence par choisir la chaine associée à la date t0 (la structure markovGenerator possède 2 générateurs, 1 pour les weekend, 1 pour les jours de semaine) (la fonction chose est dans scenarios/utils.jl)
# On chosit une date t0, on recupère les état du mois, puis on séléctionne ceux de l'heure suivant t0 et on crée un vecteur s0 - l'état pour chaque état (sur les 3 composante)
# On a alors un vecteur de différence d'état Δ_s0. Parmis lequel on cherche la norme minimum dont on récupère l'index idx_0.
# On peut résumer en : on cherche l'id de l'état du mois et de l'heure de t0 le plus proche de s0

# on stock dans t les date (timestamps) de toute les étapes de la chaines
# On peut ainsi récupérer dans hour et month, sous forme de veteur, les heures et les mois associés à ces timestamp

#On initialise les tableaux puis on passe à la génération

#On a un triple for sur nstep, ny, ns (c'est l'ordre de priorité des itérations pour un multiple for sur 1 seule ligne)
    # On séléctionne la chaine associée à l'heure h in 1:nstep
    # On selectionne pour l'état actuel (idx_0) à la "date" h la distribution de proba associé à la transition à venir
    # Le successeur est séléctionné suivant cette distribution et on stock l'état oassocié
    # idx_0 est remplacé pour faire de l'état suivant le nouvel état courat

#Pour finir on renvoi les scénario ainsi généré ainsi qu'une probabilité dont on ne se sert pas pour le moment
# Celle-ci est sensé associé à chaque scénario généré une probabilité d'existence relatve vis a vis des autres scénario
# Cette dernière probabilité est utilisé principalement pour une méthode de gestion (OLFC) se basant des scénarios à court terme pour optimiser ses opérations.
######################################################################################################################################################################################

function generate(generator::MarkovGenerator, s0, t0::DateTime, nstep::Int64; ny::Int64=1, ns::Int64=1)
    # Initialize state with the closest value
    mc = chose(generator, t0)
    Δ_s0 = [s0 .- mc.states[Dates.month(t0)][Dates.hour(t0)+1, :, :][:,state] for state in 1:generator.nstate]
    idx_0 = findmin(norm.(Δ_s0))[2]

    # Timestamp
    t = t0:Hour(1):t0+Hour(nstep-1)
    hour, month = Dates.hour.(t) .+ 1, Dates.month.(t) # hour(t[1:24]) = 0:23...

    # Pre-allocate
    ω_generated = [zeros(nstep, ny, ns) for _ in 1:size(mc.states[1],2)]
    probabilities = ones(nstep, ny, ns)

    # Initialization
    for j in 1:size(mc.states[1],2)
        ω_generated[j][1,:,:] .= s0[j]
    end

    # Generate scenarios
    for s in 1:ns, y in 1:ny, h in 1:nstep-1
        # Test to chose the appropriate markov chain
        mc = chose(generator, t[h])
        # Build categorical distribution with transition matrices - we make the assumption that the probabilities from 24h to 1h are equal whatever day transitions
        distribution = Categorical(mc.transition_matrices[month[h]][idx_0, :, hour[h]])
        # Sample state index from the distribution
        idx_1 = rand(distribution)
        # Retrieve the associated values
        for j in 1:size(mc.states[1],2)
            ω_generated[j][h+1,y,s] = mc.states[month[h]][hour[h+1], j, idx_1]
        end
        # Store probabilities
        probabilities[h+1,y,s] = probs(distribution)[idx_1]
        # Update idx
        idx_0 = idx_1
    end

    return ω_generated,   (prod(probabilities, dims=1)[1,:,:] / sum(prod(probabilities, dims=1)[1,:,:]))
end

# Perfect foresight generator
mutable struct AnticipativeGenerator <: AbstractScenariosGenerator
    forecast
    AnticipativeGenerator() = new()
end

function initialize_generator!(generator::AnticipativeGenerator, data...)
    # Store anticipative forecast
    generator.forecast = [d for d in data]
    return generator
end

# Generate perfect forecast
function generate(generator::AnticipativeGenerator, s0, t0::DateTime, nstep::Int64; ny::Int64=1, ns::Int64=1, h::Int64=1)
    # Current index
    if length(generator.forecast[1].t[h:end]) < nstep
        # Windows
        windows = h : h + length(generator.forecast[1].t[h:end]) - 1
        # Add zeros to have a constant size
        n_zeros = nstep - length(windows)
        # Forecast
        forecast = [vcat(d.power[windows], zeros(n_zeros)) for d in generator.forecast]
    else
        # Windows
        windows = h : h + nstep - 1
        # Forecast
        forecast = [d.power[windows] for d in generator.forecast]
    end

    return forecast, [1.]
end
