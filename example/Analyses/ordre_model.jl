include("..\\..\\src\\Genesys.jl")

using Main.Genesys
using JLD, Dates, Seaborn
using CSV, DataFrames
using Statistics
using StatsBase

pygui(true)

# Parameters of the simulation
const nh, ny, ns = 8760, 20, 200

# Load input data
data = load(joinpath("example","data","ausgrid_5_twostage.jld"))

#https://www.researchgate.net/figure/Battery-cycle-to-failure-versus-the-DOD-gap-34_fig5_276136726  for data1
# Modeling of Lithium-Ion Battery Degradation for Cell Life Assessment, Bolun Xu et al. for data 2
fatigue_data = DataFramesDataFrame(CSV.File("example\\data\\fatigue_data.csv.csv", delim = ";", header = [Symbol("DoD"),Symbol("cycle")], types=Dict(:DoD=>Float64, :cycle=>Float64)))


mg = Microgrid(parameters = GlobalParameters(nh, ny, ns, renewable_share = 1.))
add!(mg, Demand(carrier = Electricity()),
                Solar(),
                Liion_energy_exchanged(),
                Grid(carrier = Electricity()))

ω_d = Scenarios(mg, data["ω_optim"], adjust_length=true)



function run_simu(param::Vector{Float64}; supp_args = supp_args("vermeer", "vermeer", 2))
    controller = RBC(options = RBCOptions(policy_selection = supp_args.policy))
    designer = Manual(generations = [param[1]], storages = [param[2]], subscribed_power = [param[3]])

    couple = (E=false, R=false)

    if supp_args.soh == "linear"
        battery = Liion_energy_exchanged(;nCycle = 15000., soc_model = supp_args.soc, couplage = couple)
    elseif supp_args.soh == "rainflow"
        battery = Liion_rainflow(update_by_year = 12, soc_model = supp_args.soc, fatigue_data = fatigue_data, couplage = couple)
    elseif supp_args.soh == "electro_chimical"
        battery = Liion_electro_chimique(update_by_year = 12, soc_model = supp_args.soc, couplage = couple)
    elseif supp_args.soh == "vermeer"
        battery = Liion_vermeer(soc_model = supp_args.soc)
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

    return sum(metrics.cost.total, dims=1),mean(metrics.renewable_share[2:ny,:], dims = 1)
end


struct supp_args
    soc::String
    soh::String
    policy::Int64
end

config = [supp_args("linear", "linear", 2), supp_args("tremblay_dessaint", "linear", 2), supp_args("linear", "rainflow", 2), supp_args("tremblay_dessaint", "rainflow", 2), supp_args("linear", "electro_chimical", 2), supp_args("tremblay_dessaint", "electro_chimical", 2) ]

n_conf = length(config)

RES_matrix = zeros(n_conf, ns)
Cost_matrix = zeros(n_conf, ns)

Res_order = zeros(n_conf, ns)
Cost_order = zeros(n_conf, ns)


param = [30., 30., 10.]

for (i,conf) in enumerate(config)

    a,b = run_simu([param[1],param[2],param[3]]; supp_args = conf)
    Cost_matrix[i,:] = a[1,:]
    RES_matrix[i,:] = b[1,:]

    Cost_order[i,:] = ordinalrank(a[1,:]; lt=isless, by=identity)
    Res_order[i,:] =  ordinalrank(b[1,:]; lt=isless, by=identity)
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


old_person_cost
old_pearson_res
