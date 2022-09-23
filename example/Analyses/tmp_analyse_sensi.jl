
include("..\\..\\src\\Genesys.jl")

using Main.Genesys
using JLD, Dates, Seaborn
using CSV, DataFrames
using Statistics

pygui(true)

# Parameters of the simulation
const nh, ny, ns = 8760, 20, 1

# Load input data
data = load(joinpath("example","data","ausgrid_5_twostage.jld"))

#https://www.researchgate.net/figure/Battery-cycle-to-failure-versus-the-DOD-gap-34_fig5_276136726  for data1
# Modeling of Lithium-Ion Battery Degradation for Cell Life Assessment, Bolun Xu et al. for data 2
fatigue_data = DataFrame(CSV.File("example\\data\\fatigue_data.csv.csv", delim = ";", header = [Symbol("DoD"),Symbol("cycle")], types=Dict(:DoD=>Float64, :cycle=>Float64)))


mg = Microgrid(parameters = GlobalParameters(nh, ny, ns, renewable_share = 1.))
add!(mg, Demand(carrier = Electricity()),
                Solar(),
                Liion_energy_exchanged(),
                Grid(carrier = Electricity()))

ω_d, ω_a = Scenarios(mg, data["ω_optim"], adjust_length=true), Scenarios(mg, data["ω_simu"], adjust_length=true)


#TODO Discuter des metrics


using Sobol

mutable struct Sobol_indices
    S1::Matrix{Float64}
    ST::Matrix{Float64}
    convergence::Float64
    ci::Any #confidence_intervals
end




function generate_next_sample(s::ScaledSobolSeq, lb::Vector{Float64}, ub::Vector{Float64} )
    x = rand(Float64, length(lb)) .* (ub-lb) + lb
    return next!(s), x
end




function generate_next_observations(f, sample1, sample2, u, u_comp; supp_args = nothing)
    Y_tmp = zeros(size(sample1, 1),Nmetrics)
    Yu_tmp = zeros(size(sample1, 1), 2 * size(sample1, 2), Nmetrics)

    #kwargs = (soh = soh, soc = soc)

    for i in 1:size(sample1, 1)
        Y_tmp[i,:] = f(sample1[i,:]; supp_args = supp_args)

        for u_id in 1:(2 * size(sample1, 2))
            Yu_tmp[i,u_id,:] = f([u[u_id,1] * sample1[i,1] + u_comp[u_id,1] * sample2[i,1],
            u[u_id,2] * sample1[i,2] + u_comp[u_id,2] * sample2[i,2],
            u[u_id,3] * sample1[i,3] + u_comp[u_id,3] * sample2[i,3]]; supp_args = supp_args)
        end
    end

    return Y_tmp, Yu_tmp
end

function compute_indices(Y, Yu,  estimator::Int64)
    Nmetrics = size(Y, 2)
    Nobs = size(Y, 1)
    k = size(Yu,2)

    Z = zeros(Nobs,Nmetrics)
    M = zeros(Nobs,Nmetrics)

    #The 2 estimators
    sobol_first_id = zeros(Float64, Nvar, Nmetrics)
    sobol_tot_id =  zeros(Float64, Nvar, Nmetrics)

    for metrics in 1:Nmetrics


        #Compute indices
        for u_id in 1:k

            if estimator == 1
                numer1 = (sum(Y[i,metrics] * Yu[i,u_id,metrics] for i in 1:Nobs)/Nobs)     -    ( ( sum(Y[i,metrics] for i in 1:Nobs)/Nobs ) * ( sum(Yu[i,u_id,metrics] for i in 1:Nobs)/Nobs ))
                denom1 = (sum(Y[i,metrics]^2 for i in 1:Nobs)/Nobs)  -  ((sum(Y[i,metrics] for i in 1:Nobs)/Nobs)^2)
                if u_id <= Nvar
                    sobol_first_id[u_id,metrics] = numer1/denom1
                else
                    sobol_tot_id[u_id-Nvar,metrics] = 1 - (numer1/denom1)
                end
            elseif estimator == 2

                #Compute comon values for every index
                for i in 1:Nobs
                    Z[i,metrics] = (Y[i,metrics] + sum(Yu[i,u_id,metrics] for u_id in 1:k)) / (k+1)

                    M[i,metrics] = ( (Y[i,metrics]^2) + sum((Yu[i,u_id,metrics]^2) for u_id in 1:k)) / (k+1)
                end

                numer2 = (sum(Y[i,metrics] * Yu[i,u_id,metrics] for i in 1:Nobs)/Nobs)     -    ( sum(Y[i,metrics] + Yu[i,u_id,metrics] for i in 1:Nobs)/2Nobs )^2
                denom2 = (sum(M[i,metrics] for i in 1:Nobs)/Nobs)  -  ((sum(Z[i,metrics] for i in 1:Nobs)/Nobs)^2)

                if u_id <= Nvar
                    sobol_first_id[u_id,metrics] = numer2/denom2
                else
                    sobol_tot_id[u_id-Nvar,metrics] = 1 - (numer2/denom2)
                end
            end
        end
    end

    return sobol_first_id, sobol_tot_id
end


#Return true if indices validates conditions
function evaluate_indices(indices_over_run, indices_tot_over_run, threshold_alpha::Float64, lambda::Int64)

    threshold = 0.01
    #TODO comprendre les propriétés pour être en mesure d'écrire les tests

    #Positivity test
    if sum(indices_over_run[:,:,end] .< (0-threshold)) != 0
        return false
    end

    #positivity test
    if sum(indices_tot_over_run[:,:,end] .< (0-threshold)) != 0
        return false
    end

    #under 1 test
    if sum(indices_over_run[:,:,end]  .> (1+threshold)) != 0
        return false
    end

    #sum of first index < 1 test
    for i in 1:size(indices_over_run[:,:,end] ,2)
        if sum(indices_over_run[:,i,end]) > (1+threshold)
            return false
        end
    end

    #convergence test
    if sum(abs.(indices_over_run[:,:,(end-lambda):end] .- indices_over_run[:,:,end-(lambda+1)]) .> threshold_alpha) != 0
        return false
    end

    #All test passed
    return true

end


#Generate both necessary samples and permutations
function prepare_sampling(lb, ub, Nvar; Nobs_ini = 100)
    seq = SobolSeq(lb, ub)
    premier_sample = zeros(Float64, (Nobs_ini, Nvar))
    second_sample = zeros(Float64, (Nobs_ini, Nvar))

    #Generating the initial samples
    for row in 1:Nobs_ini
        premier_sample[row,:], second_sample[row,:] = generate_next_sample(seq, lb, ub)
    end

    #Define the combination for pick freeze method
    u_first = zeros(Int64, Nvar, Nvar)
    u_tot = ones(Int64, Nvar, Nvar)
    for i in 1:Nvar
        u_first[i,i] = 1
        u_tot[i,i] = 0
    end

    #boolean matrices used to activate the some selected variables within the the samples
    u = vcat(u_first, u_tot)
    u_comp =  ones(Int64, 2 * Nvar, Nvar) - u

    return u, u_comp, premier_sample, second_sample, seq
end



#Use a quasi monte carlo method  (with sobol seq sampling ) with variance estimators and pick freeze methode to estimates sobol indices
function get_sobol_indices(f, lb, ub ; Nobs_ini = 100, Nobs_step = 30, time_limit = 0, supp_args::Any = nothing)

    #Set the defined time limit if there's one
    # time_limit != 0 ? limit = time_from_now(time_limit) : nothing

    u, u_comp, premier_sample, second_sample, seq = prepare_sampling(lb, ub, length(ub))

    Nobs_extra = 0

    #Y is generated from the first sample, Yu contains the variations generated using the u boolean sample variable selector
    Y, Yu = generate_next_observations(f, premier_sample, second_sample, u, u_comp; supp_args = supp_args)
    S1, ST = compute_indices(Y, Yu, 2)

    #Those are computed to observe the convergence of the estimators
    indices_over_run = zeros(Float64, Nvar, 2, size(Y,1) )
    indices_tot_over_run = zeros(Float64, Nvar, 2, size(Y,1) )
    for i in 1:size(Y,1)
        indices_over_run[:,:,i], indices_tot_over_run[:,:,i] = compute_indices(Y[1:i,:], Yu[1:i,:,:], 2)
    end

    #Generate new observations while they dont fit some criterions, Nobs_step new observations are generated by loops
    while ! evaluate_indices(indices_over_run, indices_tot_over_run, 0.01, Nobs_step)


        current_Nobs = Nobs_ini + Nobs_extra

        println("observation count : ", current_Nobs)

        #extend the samples
        for row in 1:Nobs_step
            new_row1, new_row2 = generate_next_sample(seq, lb, ub)
            premier_sample = vcat(premier_sample, transpose(new_row1))
            second_sample = vcat(second_sample, transpose(new_row2))
            Nobs_extra += 1
        end

        #generate new observations with those samples
        Y_tmp, Yu_tmp = generate_next_observations(f, premier_sample[(current_Nobs+1):end,:], second_sample[(current_Nobs+1):end,:], u, u_comp; supp_args)

        #cat new observations
        Y = cat(Y,Y_tmp; dims=1)
        Yu = cat(Yu,Yu_tmp; dims=1)


        S1, ST = compute_indices(Y, Yu, 2)

        for i in (current_Nobs+1):size(Y,1)
            first, tot = compute_indices(Y[1:i,:], Yu[1:i,:,:], 2)
            indices_over_run = cat(indices_over_run, first, dims = 3)
            indices_tot_over_run = cat(indices_tot_over_run, tot, dims = 3)
        end

        #Stop if time_limit is reached
        # if time_limit != 0 && _done(limit)
        #     convergence = mean(abs.(indices_over_run[:,:,(end-Nobs_step):end] .- indices_over_run[:,:,end-(Nobs_step+1)]))
        #     return Sobol_indices(S1, ST, convergence)
        # end

    end

    convergence = mean(abs.(indices_over_run[:,:,(end-Nobs_step):end] .- indices_over_run[:,:,end-(Nobs_step+1)]))

    ci = nothing    #TODO compute confidence intervals

    return Sobol_indices(S1, ST, convergence, ci)
#   return S1, ST, Y, Yu, Nobs_extra
end

#function to determine if a time limit is reached
_done(limit) = time_ns() > limit
#funtion to get a nanosecond count with a number of minutes as input
time_from_now(minutes) = round(Int, 10^9 * 60 * minutes + time_ns())





function run_simu(param::Vector{Float64}; supp_args = supp_args("tremblay_dessaint", "rainflow", 2))
    controller = RBC(options = RBCOptions(policy_selection = supp_args.policy))
    designer = Manual(generations = [param[1]], storages = [param[2]], subscribed_power = [param[3]])

    couple = (E=false, R=false)

    if supp_args.soh == "linear"
        battery = Liion_energy_exchanged(;nCycle = 15000., soc_model = supp_args.soc, couplage = couple)
    elseif supp_args.soh == "rainflow"
        battery = Liion_rainflow(update_by_year = 12, soc_model = supp_args.soc, fatigue_file_name = "example\\data\\fatigue_data2.csv.csv", couplage = couple)
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

    return [mean(sum(metrics.cost.total,dims=1)), mean(metrics.renewable_share[2:ny,:])]
end




struct supp_args
    soc::String
    soh::String
    policy::Int64
end


############sobol indices package
using QuasiMonteCarlo, GlobalSensitivity, OrdinaryDiffEq, WGLMakie # DiffEqSensitivity

#Define main param and functions
Nobs_ini = 100

Nvar = 3
Nmetrics = 2

PV_bounds  = (min = 10., max = 50.)
Bat_bounds = (min = 10., max = 100.)
Pow_bounds = (min = 6., max = 36.)

lb = [PV_bounds.min, Bat_bounds.min, Pow_bounds.min]
ub = [PV_bounds.max, Bat_bounds.max, Pow_bounds.max]

n_var_set =  500
bounds = [[PV_bounds.min,PV_bounds.max], [Bat_bounds.min,Bat_bounds.max], [Pow_bounds.min,Pow_bounds.max]]

sampler = SobolSample()
A,B = QuasiMonteCarlo.generate_design_matrices(n_var_set,lb,ub,sampler)

#sobol_result = gsa(run_simu,GlobalSensitivity.Sobol(), A,B)






args = supp_args("linear", "linear", 2)
@time sobol_indices_lin_lin = get_sobol_indices(run_simu, lb, ub; supp_args = args)
save("sobol_indices_lin_lin_S1.jld", "sobol_indices_lin_lin_S1", sobol_indices_lin_lin.S1)
save("sobol_indices_lin_lin_ST.jld", "sobol_indices_lin_lin_ST", sobol_indices_lin_lin.ST)


args = supp_args("tremblay_dessaint", "linear", 2)
@time sobol_indices_lin_td = get_sobol_indices(run_simu, lb, ub; supp_args = args)
save("sobol_indices_lin_td_S1.jld", "sobol_indices_lin_td_S1", sobol_indices_lin_td.S1)
save("sobol_indices_lin_td_ST.jld", "sobol_indices_lin_td_ST", sobol_indices_lin_td.ST)


args = supp_args("linear", "rainflow", 2)
@time sobol_indices_rainflow_lin = get_sobol_indices(run_simu, lb, ub; supp_args = args)
save("sobol_indices_rainflow_lin_S1.jld", "sobol_indices_rainflow_lin_S1", sobol_indices_rainflow_lin.S1)
save("sobol_indices_rainflow_lin_ST.jld", "sobol_indices_rainflow_lin_ST", sobol_indices_rainflow_lin.ST)


args = supp_args("tremblay_dessaint", "rainflow", 2)
@time sobol_indices_rainflow_td = get_sobol_indices(run_simu, lb, ub; supp_args = args)
save("sobol_indices_rainflow_td_S1.jld", "sobol_indices_rainflow_td_S1", sobol_indices_rainflow_td.S1)
save("sobol_indices_rainflow_td_ST.jld", "sobol_indices_rainflow_td_ST", sobol_indices_rainflow_td.ST)


args = supp_args("linear", "electro_chimical", 2)
@time sobol_indices_elechim_lin = get_sobol_indices(run_simu, lb, ub; supp_args = args)
save("sobol_indices_elechim_lin_S1.jld", "sobol_indices_elechim_lin_S1", sobol_indices_elechim_lin.S1)
save("sobol_indices_elechim_lin_ST.jld", "sobol_indices_elechim_linST", sobol_indices_elechim_lin.ST)


args = supp_args("tremblay_dessaint", "electro_chimical", 2)
@time sobol_indices_elechim_td = get_sobol_indices(run_simu, lb, ub; supp_args = args)
save("sobol_indices_elechim_td_S1.jld", "sobol_indices_elechim_td_S1", sobol_indices_elechim_td.S1)
save("sobol_indices_elechim_td_ST.jld", "sobol_indices_elechim_td_ST", sobol_indices_elechim_td.ST)



using Plots


indices_over_run = zeros(Float64, Nvar, 2, size(Y,1) )
indices_tot_over_run = zeros(Float64, Nvar, 2, size(Y,1) )
for i in 1:size(Y,1)
    indices_over_run[:,:,i], indices_tot_over_run[:,:,i] = compute_indices(Y[1:i,:], Yu[1:i,:,:], 2)
end
Plots.plot(1:size(Y,1), transpose(indices_over_run[:,1,:]), title = "convergence of cost first index", label = ["PV" "Bat" "sous"], yaxis = [0,1] )
Plots.plot(1:size(Y,1), transpose(indices_over_run[:,2,:]),  title = "convergence of RES first index", label = ["PV" "Bat" "sous"], yaxis = [0,1] )

Plots.plot(1:size(Y,1), transpose(indices_tot_over_run[:,1,:]), title = "convergence of cost tot index", label = ["PV" "Bat" "sous"], yaxis = [0,1] )
Plots.plot(1:size(Y,1), transpose(indices_tot_over_run[:,2,:]), title = "convergence of RES tot index", label = ["PV" "Bat" "sous"], yaxis = [0,1] )


Plots.plot([maximum(abs.(indices_over_run[:,:,(i-30):i] .- indices_over_run[:,:,i-(30)])) for i in 31:size(Y,1)], yaxis=:log)


Plots.plot(transpose(sum(indices_over_run[:,1,:], dims = 1)))
Plots.plot(transpose(sum(indices_over_run[:,2,:], dims = 1)))
Plots.plot(transpose(sum(indices_tot_over_run[:,1,:], dims = 1)))
Plots.plot(transpose(sum(indices_tot_over_run[:,2,:], dims = 1)))




############THIS seems to work
using QuasiMonteCarlo, GlobalSensitivity, OrdinaryDiffEq, WGLMakie # DiffEqSensitivity


n_var_set =  Nobs_ini + Nobs_extra
bounds = [[PV_bounds.min,PV_bounds.max], [Bat_bounds.min,Bat_bounds.max], [Pow_bounds.min,Pow_bounds.max]]

sobol_result = gsa(f, GlobalSensitivity.Sobol(), bounds, N = n_var_set)









premier_sample = zero(Float64, (Nobs, Nvar))
second_sample = rand(Float64, (Nobs, Nvar))

premier_sample[:,1] = premier_sample[:,1] * (PV_bounds.max - PV_bounds.min) .+ PV_bounds.min
premier_sample[:,2] = premier_sample[:,2] * (Bat_bounds.max - Bat_bounds.min) .+ Bat_bounds.min
premier_sample[:,3] = premier_sample[:,3] * (Pow_bounds.max - Pow_bounds.min) .+ Pow_bounds.min

second_sample[:,1] = second_sample[:,1] * (PV_bounds.max - PV_bounds.min) .+ PV_bounds.min
second_sample[:,2] = second_sample[:,2] * (Bat_bounds.max - Bat_bounds.min) .+ Bat_bounds.min
second_sample[:,3] = second_sample[:,3] * (Pow_bounds.max - Pow_bounds.min) .+ Pow_bounds.min








mutable struct Labels_battery
     policy::String
     SOC::String
     SOH::String
     couplageE::Bool
     couplageR::Bool
end






#List of explored combinations
SOH_list = ["linear", "rainflow", "electro_chimical"]
SOC_list = ["linear", "tremblay_dessaint"]
couplage = [(E = false, R = false), (E = true, R = false),(E = false, R = true),(E = true, R = true)]




############THIS seems to work


##########################
# Sensistibity analysis #
##########################
     labels_list = []
     n_var_set = 1000
     bounds = [[PV_bounds.min,PV_bounds.max], [Bat_bounds.min,Bat_bounds.max], [Pow_bounds.min,Pow_bounds.max]]
     n_config = length(SOC_list) * length(SOH_list) * length(couplage)


     cost_Total_order_indices = Matrix{Any}(undef, n_config, 3)
     cost_first_order_indices = Matrix{Any}(undef, n_config, 3)
     RES_Total_order_indices = Matrix{Any}(undef, n_config, 3)
     RES_first_order_indices = Matrix{Any}(undef, n_config, 3)

     i_config = 1

     for soh in SOH_list
         for soc in SOC_list
            # for couple in couplage

                 #apply battery params
                 if soh == "linear"
                     battery = Liion_energy_exchanged(;nCycle = 15000., soc_model = soc, couplage = couple)
                 elseif soh == "rainflow"
                     battery = Liion_rainflow(update_by_year = 12, soc_model = soc, fatigue_file_name = "example\\data\\fatigue_data2.csv.csv", couplage = couple)
                 elseif soh == "electro_chimical"
                     battery = Liion_electro_chimique(update_by_year = 12, soc_model = soc, couplage = couple)
                 end

                 label = Labels_battery(string(2), soh, soc, couple.E, couple.R)
                 @time begin
                     sobol_result = gsa(f, Sobol(), bounds, N = n_var_set)
                 end
                 cost_Total_order_indices[i_config,:] = sobol_result.ST[2,:]
                 cost_first_order_indices[i_config,:] = sobol_result.S1[2,:]
                 RES_Total_order_indices[i_config,:] = sobol_result.ST[1,:]
                 RES_first_order_indices[i_config,:] = sobol_result.S1[1,:]

                 push!(labels_list, label)

                 i_config += 1

             #end
        end
    end
gsa(f,eFAST(),bounds,n=n_var_set)


pygui(true)

plot( (data["ω_optim"]["ld_E"].power[:,1,1]+data["ω_optim"]["ld_H"].power[:,1,1]) .- data["ω_optim"]["pv"].power[:,1,1])

plot( data["ω_optim"]["ld_E"].power[:,1,:] .- data["ω_optim"]["pv"].power[:,1,:])


std_diff = zeros(8760)
std_mean_diff = zeros(8760)

std_diff_day = zeros(365)
std_mean_diff_day = zeros(365)

for h in 1:8760
    Pdiff = (data["ω_optim"]["ld_E"].power[h,y,:]+data["ω_optim"]["ld_H"].power[h,y,:]) .- data["ω_optim"]["pv"].power[h,y,:]
    std_diff[h] = std(Pdiff)
    std_mean_diff[h] = std(Pdiff) / mean(Pdiff)
end

for d in 1:(365)
    h = ((d-1) * 24 + 1):(d * 24)
    Pdiff = sum((data["ω_optim"]["ld_E"].power[h,y,:]+data["ω_optim"]["ld_H"].power[h,y,:]) .- data["ω_optim"]["pv"].power[h,y,:], dims = 1)
    std_diff_day[d] = std(Pdiff)
    std_mean_diff_day[d] = std(Pdiff) / mean(Pdiff)
end


plot(std_diff)
plot(std_mean_diff)

plot(std_diff_day)
plot(std_mean_diff_day)
