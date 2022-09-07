# Load packages
include("..\\src\\Genesys.jl")

using Main.Genesys
using JLD, Dates#, Seaborn
using CSV, DataFrames

using Statistics
using GaussianProcesses
using Random
using Plots
#using PyPlot



pygui(true)


# Parameters of the simulation
const nh, ny, ns = 8760, 20, 1

# Load input data
data = load(joinpath("example","data","ausgrid_5_twostage.jld"))

microgrid_non_linear_asset = Microgrid(parameters = GlobalParameters(nh, ny, ns, renewable_share = 0.9))

microgrid_linear_asset = Microgrid(parameters = GlobalParameters(nh, ny, ns, renewable_share = 0.9))

add!(microgrid_non_linear_asset, Demand(carrier = Electricity()),
                Solar(),
                Liion_electro_chimique(update_by_year = 12, soc_model = "tremblay_dessaint"),
                Grid(carrier = Electricity()))


add!(microgrid_linear_asset, Demand(carrier = Electricity()),
                Solar(),
                Liion_energy_exchanged(),
                Grid(carrier = Electricity()))


                # Initialize scenarios

ω_d = Scenarios(microgrid_non_linear_asset, data["ω_optim"], rep = 10)
#ω_a = Scenarios(microgrid_simulation, data["ω_optim"], rep = 10), Scenarios(microgrid_simulation, data["ω_simu"], rep=10)


costs = []
share = []
opex = []
capex = []
salvage = []

N = 2000

max_PV = 100
min_PV = 5
max_B = 100
min_B = 5
x = rand(2, N)
x[1,:] = (x[1,:] * (max_PV-min_PV)) .+ min_PV #Between min_PV and max_PV
x[2,:] = (x[2,:] * (max_B-min_B)) .+ min_B    #Between min_B and max_B

model_share = 0.75

feed_model_interval = 1:convert(Int64,(model_share * N))
test_model_interval = convert(Int64,(model_share * N)+1):N

#Have a look to the coverage
df = DataFrame(Batterie=x[1,:], PV=x[2,:], Type = vcat(repeat(["feed"],length(feed_model_interval)), repeat(["test"],length(test_model_interval))))
using PlotlyJS
PlotlyJS.plot(
    df, x=:Batterie, y=:PV, color=:Type,
    marker=attr(size=5, sizemode="area"),
    mode="markers"
)

#Simulate N time to get N observations
for i in 1:N
    designer_non_linear = initialize_designer!(microgrid_non_linear_asset, Manual(generations = [x[1,i]], storages = [x[2,i]]), ω_d) #Manual(generations = [x[1,i]], storages = [x[2,i]]), ω_d)

    controller_non_linear = initialize_controller!(microgrid_non_linear_asset, RBC(options = RBCOptions(policy_selection = 2)), ω_d)

    simulate!(microgrid_non_linear_asset, controller_non_linear, designer_non_linear, ω_d, options = Genesys.Options(mode = "serial", firstyear = false))
    metrics_non_linear = Metrics(microgrid_non_linear_asset, designer_non_linear)

    push!(costs, sum(metrics_non_linear.npv.capex[:,1] + metrics_non_linear.npv.opex[:,1] - metrics_non_linear.npv.salvage[:,1]))
    push!(share, mean(metrics_non_linear.renewable_share[2:ny,1]))
    push!(opex,  sum(metrics_non_linear.npv.opex[:,1]))
    push!(capex, sum(metrics_non_linear.npv.capex[:,1]))
    push!(salvage, sum(metrics_non_linear.npv.salvage[:,1]))


    println(i,"/",N)
end



MAE_cost = []
MAE_share = []
MSE_cost = []
MSE_share = []
RE_cost = []
RE_share = []
MinE_cost = []
MinE_share = []
MaxE_cost = []
MaxE_share = []


Mean_costs_list = []
Mean_share_list = []

push!(Mean_costs_list, MeanZero())
push!(Mean_share_list, MeanZero())

Kern_cost_list = []
Kern_share_list = []

push!(Kern_cost_list, SE(0.5,0.0))
push!(Kern_cost_list, SE(0.5,0.0) * Lin(0.0))
push!(Kern_cost_list, Matern(3/2,0.0,0.0))


push!(Kern_share_list, SE(0.5,0.0))
push!(Kern_share_list, SE(0.5,0.0) * Lin(0.0))
push!(Kern_share_list, Matern(3/2,0.0,0.0))




logObsNoise = -2.0                        # log standard deviation of observation noise (this is optional)


using Optim

for i in 1:length(Mean_costs_list)
    for j in 1:length(Kern_cost_list)
        println("i = ", i, " j = ", j)
        gp_costs = GP(x[:,feed_model_interval], costs[feed_model_interval], Mean_costs_list[i], Kern_cost_list[j], logObsNoise)
        gp_share = GP(x[:,feed_model_interval], share[feed_model_interval], Mean_share_list[i], Kern_share_list[j], logObsNoise)

        optimize!(gp_costs)
        optimize!(gp_share)

        μ, σ² = predict_y(gp_costs, x[:, test_model_interval])

        push!(MAE_cost, (1/length(test_model_interval)) * sum(abs.(μ - costs[test_model_interval])))
        push!(MSE_cost, (1/length(test_model_interval)) * sum(abs.(μ - costs[test_model_interval]).^2))
        push!(RE_cost, (1/length(test_model_interval)) * sum(abs.((μ - costs[test_model_interval])./costs[test_model_interval])) *100)
        push!(MinE_cost, minimum(abs.((μ - costs[test_model_interval])./costs[test_model_interval]))*100)
        push!(MaxE_cost,  maximum(abs.((μ - costs[test_model_interval])./costs[test_model_interval]))*100)

        μ, σ² = predict_y(gp_share,x[:, test_model_interval])

        push!(MAE_share, (1/length(test_model_interval)) * sum(abs.(μ - share[test_model_interval])))
        push!(MSE_share, (1/length(test_model_interval)) * sum(abs.(μ - share[test_model_interval]).^2))
        push!(RE_share, (1/length(test_model_interval)) * sum(abs.((μ - share[test_model_interval])./share[test_model_interval]))*100)
        push!(MinE_share, minimum(abs.((μ - share[test_model_interval])./share[test_model_interval]))*100)
        push!(MaxE_share,  maximum(abs.((μ - share[test_model_interval])./share[test_model_interval]))*100)

    end
end

selected_model_cost = sortperm(RE_cost)[1]
selected_model_share = sortperm(RE_share)[1]
gp_costs = GP(x[:,feed_model_interval], costs[feed_model_interval], Mean_costs_list[1], Kern_cost_list[selected_model_cost], logObsNoise)
optimize!(gp_costs)
Plots.surface(gp_costs, title = "cost", xlabel="PV", ylabel = "battery")
Plots.heatmap(gp_costs, title = "cost", xlabel="PV", ylabel = "battery")


gp_share = GP(x[:,feed_model_interval], share[feed_model_interval], Mean_share_list[1], Kern_share_list[selected_model_share], logObsNoise)
optimize!(gp_share)
Plots.surface(gp_share, title = "share", xlabel="PV", ylabel = "battery")
Plots.heatmap(gp_share, title = "share", xlabel="PV", ylabel = "battery")
