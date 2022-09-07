include("..\\src\\Genesys.jl")

using Main.Genesys
using Distributions, CSV, JLD, Dates, Seaborn, Statistics, ProgressMeter

#Le package pour les algo de clustering
using Clustering 

pygui(true)

# Global constant
const nh, ny, ns = 8760, 1, 2000

#La structure qu'on utilise pour stocker les etats et transitions
mutable struct MarkovChain
    states
    transition_matrices
end


#Le format de date
dateformat = Dates.DateFormat("dd-u-yyyy")

#On charge les données relative à plusieurs clients sur 3 ans.
data_2010 = load(joinpath("stage_scenario", "data", "clean_dataset_ausgrid_2010_2011.jld"), "clean_dataset")
data_2011 = load(joinpath("stage_scenario", "data", "clean_dataset_ausgrid_2011_2012.jld"), "clean_dataset")
data_2012 = load(joinpath("stage_scenario", "data", "clean_dataset_ausgrid_2012_2013.jld"), "clean_dataset")

#On trouve l'ID des clients qui ont des données sur ces 3 ans
# Find custumer with 3 years of data
custumers = [k for k in keys(data_2010) if k in keys(data_2011) && k in keys(data_2012)]

#Je séléctionne un client
c = custumers[1]

#Je concatène les données sur 3 ans pour chaque type de données
#éléctrique
ld_E = (t = hcat(Dates.DateTime.(data_2010[c]["GC"]["time"][1:nh], dateformat), Dates.DateTime.(data_2011[c]["GC"]["time"][1:nh], dateformat), Dates.DateTime.(data_2012[c]["GC"]["time"][1:nh], dateformat)),
        power = hcat(data_2010[c]["GC"]["power"][1:nh], data_2011[c]["GC"]["power"][1:nh], data_2012[c]["GC"]["power"][1:nh]))
#heat
ld_H = (t = hcat(Dates.DateTime.(data_2010[c]["CL"]["time"][1:nh], dateformat), Dates.DateTime.(data_2011[c]["CL"]["time"][1:nh], dateformat), Dates.DateTime.(data_2012[c]["CL"]["time"][1:nh], dateformat)),
        power = hcat(data_2010[c]["CL"]["power"][1:nh], data_2011[c]["CL"]["power"][1:nh], data_2012[c]["CL"]["power"][1:nh]))
#PV
pv = (t = hcat(Dates.DateTime.(data_2010[c]["GG"]["time"][1:nh], dateformat), Dates.DateTime.(data_2011[c]["GG"]["time"][1:nh], dateformat), Dates.DateTime.(data_2012[c]["GG"]["time"][1:nh], dateformat)),
        power = hcat(data_2010[c]["GG"]["power"][1:nh], data_2011[c]["GG"]["power"][1:nh], data_2012[c]["GG"]["power"][1:nh]))



#J'agglomère les 3 type en un seul vecteur qui contient 3 tuple {power, t} où power et t sont des séries de Datetimes et de Float
data = [ld_E, ld_H, pv]


function clustering_month(data, t; flag::String="")
    # Parameters
    nm = 12 # number of month in one year

    # Pre_allocate
    data_cluster= Vector{Array{Float64,2}}(undef, nm)

    for m in 1:nm
        data_cluster[m] = reshape(data[(Dates.month.(t) .== m)], 24, :)
    end

    return data_cluster
end



# Clustering par month pour chaque  "data series" (pv, ld_E, ld_H)
months_data = [clustering_month(d.power, d.t) for d in data]


#On choisi un algo de clusterisation
algo = "kmeans"
# et le nombre de cluster qui va avec
nstate = 20




nk = size(months_data,1) # number of data series (3)
nh = 24 # number of hours (24)
nm = size(months_data[1],1) # number of month (12)

# Pre-allocate (on prépare les matrices aceuillant nos états (nh*nk*nstate) pour chaque mois )
states = [zeros(nh, nk, nstate) for m in 1:nm]
#On prépqre les matrices de transitions pour chaque mois  (nstate*nstate* nh-1)
# nh-1 car on ne déféni pas la transition de la dernière à la première heure. C'est d'ailleurs l'un des points que l'on souhaite étudier
transition_matrices = [zeros(nstate, nstate, nh - 1) for m in 1:nm]

# Pour chaque mois on rempli les matrices d'état
for m in 1:nm
    data_for_1_month = [months_data[k][m] for k in 1:nk]

    # Normalization
    data_n = replace!.(data_for_1_month ./ maximum.(data_for_1_month), NaN => 0.)

    # number of available days for this month and type of day called N ou N_obs
    nd = size(months_data[1][m], 2)

    # On défini les enchainement de données de cluster en clsuter selon les donnéess orginales
    # Elle permettent ensuite de calculer les probabilités de transistion.
    # Ici on prédéfini simplement la matrices
    sequence = zeros(Int64, nh+1, nd)

        #########################################
        # COMPUTE STATES
        #########################################

        for h in 1:nh

            # Aggregate normalized data
            data_agg = permutedims(hcat([data_n[k][h,:] for k in 1:nk]...))
            # Clustering (ici j'ai mis le K-means mais c'est pareil pour les autres)
            clusters = kmeans(data_agg, nstate)

            # Store the denormalized states and sequences
            states[m][h,:,:] = clusters.centers .* maximum.(data_for_1_month)

            #ici on la rempli grace à l'assignement des data points dans les différents clusters
            sequence[h,:] = permutedims(clusters.assignments)

            if h == 1
                # ici on souhaite remettre en question ce choix
                sequence[end, :] = sequence[1, :]
            end
        end
        #############################################
        # COMPUTE MATRICES DE TRANSITIONS
        #############################################


        nstep = size(sequence, 1) # number of step in the chain

        # Pre-allocate
        transition_matrix = zeros(nstate, nstate, nstep - 1)

        for h in 1:nh-1
            #Ici on incrémente les transistion qui existe de cluster a cluster pour définir les proba de transitions
            for (i,j) in zip(sequence[h, :], sequence[h+1, :]) #https://docs.julialang.org/en/v1/base/iterators/
                transition_matrix[i, j, h] += 1
            end

            # On transforme ces sommes en probabilité
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

end


#On crée notre structure chaine de markov que l'on utilisera pour la génération.
# Celle-ci contient les états, ainsi que les matrices de transitions.
MC = MarkovChain(states, transition_matrices)
