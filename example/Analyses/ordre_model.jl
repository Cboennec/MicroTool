include("..\\..\\src\\Genesys.jl")

using Main.Genesys
using JLD, Dates, Seaborn
using CSV, DataFrames
using Statistics
using StatsBase

pygui(true)

# Parameters of the simulation
const nh, ny, ns = 8760, 21, 2

# Load input data
data = load(joinpath("example","data","ausgrid_5_twostage.jld"))

#https://www.researchgate.net/figure/Battery-cycle-to-failure-versus-the-DOD-gap-34_fig5_276136726  for data1
# Modeling of Lithium-Ion Battery Degradation for Cell Life Assessment, Bolun Xu et al. for data 2
fatigue_data = DataFrame(CSV.File("example\\data\\fatigue_data_NMC.csv", delim = ",", header = [Symbol("DoD"),Symbol("cycle")], types=Dict(:DoD=>Float64, :cycle=>Float64)))


mg = Microgrid(parameters = GlobalParameters(nh, ny, ns, renewable_share = 1.))
add!(mg, Demand(carrier = Electricity()),
                Solar(),
                Liion_energy_exchanged(),
                Grid(carrier = Electricity()))

ω_d = Scenarios(mg, data["ω_optim"]; same_year=false, seed=1:1000)



function run_simu(param::Vector{Float64}; supp_args = supp_args("fixed_lifetime", "fixed_lifetime", 2, true))
    controller = RBC(options = RBCOptions(policy_selection = supp_args.policy))
    designer = Manual(generations = [param[1]], storages = [param[2]], subscribed_power = [param[3]])

    couple = supp_args.couple ? (E=true, R=true) : (E=false, R=false)


    if supp_args.soh == "linear"
        battery = Liion_energy_exchanged(;calendar = true, nCycle = fatigue_data.cycle[findfirst(fatigue_data.DoD .> (0.6))], soc_model = supp_args.soc, couplage = couple)
    elseif supp_args.soh == "rainflow"
        battery = Liion_rainflow(update_by_year = 12, calendar = true, soc_model = supp_args.soc, fatigue_data = fatigue_data, couplage = couple)
    elseif supp_args.soh == "electro_chimical"
        battery = Liion_electro_chimique(update_by_year = 12, soc_model = supp_args.soc, couplage = couple)
    elseif supp_args.soh == "fixed_lifetime"
        battery = Liion_fixed_lifetime(soc_model = supp_args.soc,  couplage = couple )
    else
        error("non existing value ", soh, " for kwargs soh in function 'f'")
    end

    mg = Microgrid(parameters = GlobalParameters(nh, ny, ns, renewable_share = 1.))

    #Create microgrid
    add!(mg, Demand(carrier = Electricity()),
                    Solar(),
                    battery,
                    Grid(carrier = Electricity()))

    #generate scenarios  #Too many for the moment by the way, the need is different for this study

    designer = initialize_designer!(mg, designer, ω_d)
    controller = initialize_controller!(mg, controller, ω_d)

    simulate!(mg, controller, designer, ω_d, options = Genesys.Options(mode = "multithreads", firstyear = false))

    metrics = Metrics(mg, designer)

    println(string(typeof(mg.storages[1])))

    return sum(metrics.cost.total, dims=1),mean(metrics.renewable_share[2:ny,:], dims = 1)
end


mutable struct supp_args
    soc::String
    soh::String
    policy::Int64
    couple::Bool
end

config = [supp_args("linear", "fixed_lifetime", 2, false), supp_args("linear", "fixed_lifetime", 2, true),
  supp_args("linear", "linear", 2, false ), supp_args("linear", "linear", 2, true),
  supp_args("linear", "rainflow", 2, false), supp_args("linear", "rainflow", 2, true),
  supp_args("linear", "electro_chimical", 2, false), supp_args("linear", "electro_chimical", 2, true)]

n_conf = length(config)



using Sobol

PV_range = [5,50]
BAT_range = [10,80]

P_sous = 10.

Nb_sizing1 = 256

lb = [PV_range[1], BAT_range[1], P_sous]
ub = [PV_range[2], BAT_range[2], P_sous]
seq = SobolSeq(lb, ub)
sizings = transpose(Base.reduce(hcat, next!(seq) for i = 1:Nb_sizing1))
sizings_tuple1 = [(pv = sizings[i,1], bat = sizings[i,2], P_sub = sizings[i,3]) for i in 1:Nb_sizing1]

Seaborn.scatter(sizings[:,1], sizings[:,2])

PV_range = [5,50]
BAT_range = [10,80]
P_sous = 10.
Nb_sizing2 = 256
rnd = rand(2, Nb_sizing2)
rnd[1,:] = rnd[1,:] * (PV_range[2] - PV_range[1]) .+ PV_range[1]
rnd[2,:] = rnd[2,:] * (BAT_range[2] - BAT_range[1]) .+ BAT_range[1]

sizings = hcat(transpose(rnd), 10*ones(Nb_sizing2))
sizings_tuple2 = [(pv = sizings[i,1], bat = sizings[i,2], P_sub = sizings[i,3]) for i in 1:Nb_sizing2]
Seaborn.scatter(sizings[:,1], sizings[:,2])

sizings_tuple = cat(sizings_tuple1, sizings_tuple2, dims=1)


Nb_sizing = Nb_sizing1 + Nb_sizing2


RES_matrix = zeros(n_conf, Nb_sizing)
Cost_matrix = zeros(n_conf, Nb_sizing)

Res_order = zeros(n_conf, Nb_sizing)
Cost_order = zeros(n_conf, Nb_sizing)


for (i,conf) in enumerate(config)
    for j in 1:Nb_sizing
        a,b = run_simu([sizings_tuple[j].pv, sizings_tuple[j].bat, sizings_tuple[j].P_sub]; supp_args = conf)
        Cost_matrix[i,j] = a[1,1]
        RES_matrix[i,j] = b[1,1]

    end
    Cost_order[i,:] = ordinalrank(Cost_matrix[i,:]; lt=isless, by=identity)
    Res_order[i,:] =  ordinalrank(RES_matrix[i,:]; lt=isless, by=identity)
end



function get_pearson_coef(order_a, order_b)
    return rs = 1 - ((6 * sum( (order_a[i] - order_b[i])^2 for i in 1:length(order_b) ))/(length(order_b)*(length(order_b)^2-1)))
end


pearson_matrix_res = ones(n_conf, n_conf)
pearson_matrix_cost = ones(n_conf, n_conf)

for i in 1:length(config)
    for j in 1:length(config)

        pearson_matrix_res[i,j] = get_pearson_coef(Res_order[i,:], Res_order[j,:])
        pearson_matrix_cost[i,j] = get_pearson_coef(Cost_order[i,:], Cost_order[j,:])

    end
end


pearson_matrix_res
pearson_matrix_cost

round.(pearson_matrix_cost, digits = 2)
