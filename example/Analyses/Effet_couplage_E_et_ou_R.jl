# On souhaite ici mesurer si le couplage E et R ont un effet (individuel et collectif) sur la dégradation de la batterie et l'autonomie RES, si cela dépend exliquer de quoi cela dépend

#hypothèse : aucun impact sur les petites batterie, un impact assez important sur les grosses batteries.
            #Dans la fourchette 100-80% tout sela est encore plus négligeable



# Observer les effets de E, observer les effets de R, observer l'effet croisé
# Pour chaque jeux de donnée faire tourner les 4 configurations.

#On cherche a voir si une différence existe sur le coût général, l'autonomie, la durée de vie de la batterie.
#On test alors sur de nombreux jeux de données les 4 configs.
#On fait pour les 3 métriques , 3 heat map de la différence en % entre 0 couplage et les 3 couplages.

#Ce qui fait 3 * 3 heatmap

include("..\\..\\src\\Genesys.jl")

using Main.Genesys
using JLD, Dates, Seaborn
using CSV, DataFrames
using Statistics
using Plots, ColorSchemes


pygui(true)

# Parameters of the simulation
const nh, ny, ns = 8760, 20, 1

# Load input data
data = load(joinpath("example","data","ausgrid_5_twostage.jld"))

#https://www.researchgate.net/figure/Battery-cycle-to-failure-versus-the-DOD-gap-34_fig5_276136726  for data1
# Modeling of Lithium-Ion Battery Degradation for Cell Life Assessment, Bolun Xu et al. for data 2
fatigue_data = DataFrame(CSV.File("example\\data\\fatigue_data.csv.csv", delim = ";", header = [Symbol("DoD"),Symbol("cycle")], types=Dict(:DoD=>Float64, :cycle=>Float64)))



function get_battery_lifetime(liion::AbstractLiion)

	nb_bat = 0

	for i in 1:(ny+1)
		if liion.soh[end,i,1] <= liion.SoH_threshold
			nb_bat += 1
		end
	end

	if liion.soh[1,end,1] > liion.SoH_threshold
		nb_bat += (1-liion.soh[1,end,1])/(1 - liion.SoH_threshold)
	end

	return ny/nb_bat
end


microgrid_fixed = Microgrid(parameters = GlobalParameters(nh, ny, ns, renewable_share = 1.))
add!(microgrid_fixed, Demand(carrier = Electricity()),
                Solar(),
                Liion_fixed_lifetime(soc_model = "linear", lifetime=2),
                Grid(carrier = Electricity()))

ω_d = Scenarios(microgrid_fixed, data["ω_optim"], adjust_length=true)


param_rangepv =  5:5:20
param_rangebat = 5:5:10
param_range = 5:10
df1_result = DataFrame(pv_size = Int[], bat_size = Int[], couple= NamedTuple{(:E, :R), Tuple{Bool, Bool}}[], lifetime=Float64[], RES=Float64[], cost=Float64[])
df2_result = DataFrame(pv_size = Int[], bat_size = Int[], couple= NamedTuple{(:E, :R), Tuple{Bool, Bool}}[], lifetime=Float64[], RES=Float64[], cost=Float64[])
df3_result = DataFrame(pv_size = Int[], bat_size = Int[], couple= NamedTuple{(:E, :R), Tuple{Bool, Bool}}[], lifetime=Float64[], RES=Float64[], cost=Float64[])


for model in 1:3
	if model == 1

		it = 0
		for bat_size in param_rangebat
			for pv_size in param_rangepv
			it+=1
			println(it, "/", length(param_range)^2)

			param = [pv_size, bat_size, 5]
			controller = RBC(options = RBCOptions(policy_selection = 2 ))
			designer = Manual(generations = [param[1]], storages = [param[2]], subscribed_power = [param[3]])

				for E in [true,false]
					for R in [true, false]
						battery = Liion_energy_exchanged(;nCycle = 15000., soc_model = "tremblay_dessaint", couplage = (E=E,R=R))

						#Liion_rainflow(update_by_year = 12, fatigue_file_name = "example\\data\\fatigue_data.csv.csv" , soc_model = "tremblay_dessaint", couplage = (E=E,R=R))

						mg = Microgrid(parameters = GlobalParameters(nh, ny, ns, renewable_share = 1.))

						#Create microgrid
						add!(mg, Demand(carrier = Electricity()),
						                Solar(),
						                battery,
						                Grid(carrier = Electricity()))

						#ω_d = Scenarios(mg, data["ω_optim"], adjust_length=true)
						designer = initialize_designer!(mg, designer, ω_d)
						controller = initialize_controller!(mg, controller, ω_d)

						simulate!(mg, controller, designer, ω_d, options = Genesys.Options(mode = "serial", firstyear = false))

						lifetime = get_battery_lifetime(mg.storages[1])
						metrics = Metrics(mg, designer)
		                RES = mean(metrics.renewable_share[2:end])
						cost = sum(metrics.cost.total)
						push!(df1_result, (pv_size, bat_size, (E=E,R=R), lifetime, RES, cost))
					end
				end
			end
		end


	elseif model == 2

		it = 0
		for bat_size in param_rangebat
			for pv_size in param_rangepv
			it+=1
			println(it, "/", length(param_range)^2)

			param = [pv_size, bat_size, 5]
			controller = RBC(options = RBCOptions(policy_selection = 2 ))
			designer = Manual(generations = [param[1]], storages = [param[2]], subscribed_power = [param[3]])

				for E in [true,false]
					for R in [true, false]
						battery =   Liion_rainflow(soc_model = "tremblay_dessaint", couplage = (E=E,R=R))

						#Liion_rainflow(update_by_year = 12, fatigue_file_name = "example\\data\\fatigue_data.csv.csv" , soc_model = "tremblay_dessaint", couplage = (E=E,R=R))

						mg = Microgrid(parameters = GlobalParameters(nh, ny, ns, renewable_share = 1.))

						#Create microgrid
						add!(mg, Demand(carrier = Electricity()),
										Solar(),
										battery,
										Grid(carrier = Electricity()))

						#ω_d = Scenarios(mg, data["ω_optim"], adjust_length=true)

						designer = initialize_designer!(mg, designer, ω_d)
						controller = initialize_controller!(mg, controller, ω_d)

						simulate!(mg, controller, designer, ω_d, options = Genesys.Options(mode = "serial", firstyear = false))

						lifetime = get_battery_lifetime(mg.storages[1])
						metrics = Metrics(mg, designer)
						RES = mean(metrics.renewable_share[2:end])
						cost = sum(metrics.cost.total)
						push!(df2_result, (pv_size, bat_size, (E=E,R=R), lifetime, RES, cost))
					end
				end
			end
		end



	elseif model == 3
		it = 0
		for bat_size in param_range
			for pv_size in param_rangepv
			it+=1
			println(it, "/", length(param_range)^2)

			param = [pv_size, bat_size, 5]
			controller = RBC(options = RBCOptions(policy_selection = 2 ))
			designer = Manual(generations = [param[1]], storages = [param[2]], subscribed_power = [param[3]])

				for E in [true,false]
					for R in [true, false]
						battery =  Liion_electro_chimique(soc_model = "tremblay_dessaint", couplage = (E=E,R=R))

						#Liion_rainflow(update_by_year = 12, fatigue_file_name = "example\\data\\fatigue_data.csv.csv" , soc_model = "tremblay_dessaint", couplage = (E=E,R=R))

						mg = Microgrid(parameters = GlobalParameters(nh, ny, ns, renewable_share = 1.))

						#Create microgrid
						add!(mg, Demand(carrier = Electricity()),
										Solar(),
										battery,
										Grid(carrier = Electricity()))

						#ω_d = Scenarios(mg, data["ω_optim"], adjust_length=true)

						designer = initialize_designer!(mg, designer, ω_d)
						controller = initialize_controller!(mg, controller, ω_d)

						simulate!(mg, controller, designer, ω_d, options = Genesys.Options(mode = "serial", firstyear = false))

						lifetime = get_battery_lifetime(mg.storages[1])
						metrics = Metrics(mg, designer)
						RES = mean(metrics.renewable_share[2:end])
						cost = sum(metrics.cost.total)
						push!(df3_result, (pv_size, bat_size, (E=E,R=R), lifetime, RES, cost))
					end
				end
			end
		end

	end





end

#FIND Boundaries


val_min = ones(3,3,3).* Inf #i model, j metrique
val_max = ones(3,3,3).* -Inf

for (model, df_result) in enumerate([df1_result, df2_result, df3_result])

	df_ = filter(:couple => ==((E=false,R=false)), df_result)
	df_E = filter(:couple => ==((E=true,R=false)), df_result)
	df_R = filter(:couple => ==((E=false,R=true)), df_result)
	df_ER = filter(:couple => ==((E=true,R=true)), df_result)

	diff_lifetime = DataFrame(E = (df_E.lifetime .- df_.lifetime)./df_.lifetime, R = (df_R.lifetime .- df_.lifetime)./df_.lifetime, ER= (df_ER.lifetime .- df_.lifetime)./df_.lifetime)
	diff_RES = DataFrame(E = (df_E.RES .- df_.RES)./df_.RES, R = (df_R.RES .- df_.RES)./df_.RES, ER= (df_ER.RES .- df_.RES)./df_.RES)
	diff_cost = DataFrame(E = (df_E.cost .- df_.cost)./df_.cost, R = (df_R.cost .- df_.cost)./df_.cost, ER= (df_ER.cost .- df_.cost)./df_.cost)

	for (j, data) in enumerate([diff_cost, diff_RES, diff_lifetime])
		for (k,couplage) in enumerate([:E,:R,:ER])
			 val_min[model, j, k] = minimum(data[!,couplage])
			 val_max[model, j, k] = maximum(data[!,couplage])
			#val_min[model,k] = min(minimum(data[!,:E]),minimum(data[!,:R]),minimum(data[!,:ER]), val_min[i,j])
			#val_max[model,k] = max(maximum(data[!,:E]),maximum(data[!,:R]),maximum(data[!,:ER]), val_max[i,j])
		end
	end
end

#Plot result


mycolors = [:balance, :solar, :tofino]

for model in 1:3
	if model == 2

		folder_name = string(typeof(Liion_rainflow()))

		df_ = filter(:couple => ==((E=false,R=false)), df2_result)
		df_E = filter(:couple => ==((E=true,R=false)), df2_result)
		df_R = filter(:couple => ==((E=false,R=true)), df2_result)
		df_ER = filter(:couple => ==((E=true,R=true)), df2_result)

		diff_lifetime = DataFrame(E = (df_.lifetime .- df_E.lifetime)./df_.lifetime, R = (df_.lifetime .- df_R.lifetime)./df_.lifetime, ER= (df_.lifetime .- df_ER.lifetime)./df_.lifetime)
		diff_RES = DataFrame(E = (df_.RES .- df_E.RES)./df_.RES, R = (df_.RES .- df_R.RES)./df_.RES, ER= (df_.RES .- df_ER.RES)./df_.RES)
		diff_cost = DataFrame(E = (df_.cost .- df_E.cost)./df_.cost, R = (df_.cost .- df_R.cost)./df_.cost, ER= (df_.cost .- df_ER.cost)./df_.cost)




		x = [string(unique(df_.pv_size)[i]) for i = 1:length(unique(df_.pv_size))]
		y = [string(unique(df_.bat_size)[i]) for i = 1:length(unique(df_.pv_size))]


		if !isdir(string("figures\\", folder_name))
			mkdir(string("figures\\", folder_name))
		end







		i = 0
		for (metric, data) in enumerate([diff_cost, diff_RES, diff_lifetime])
			i += 1

			if i == 1
				abreviation = "L"
				metric_name = "Lifetime"
			elseif i == 2
				abreviation = "RES"
				metric_name = "RES"
			elseif i == 3
				abreviation = "C"
				metric_name = "Cost"

			end

			for (k,type) in enumerate([:E,:R,:ER])
				z=reshape(data[!,type], length(unique(df_.bat_size)), length(unique(df_.bat_size)))
				Plots.heatmap(x, y, abs.(z*100) , aspect_ratio = 0.7,
				 title = string(metric_name ," impact of ",type," couplage"), xlabel = "PV size", ylabel = "Battery size",
				 colorbar_title = string("\n\n\n(",abreviation, " - ", abreviation, "')/", abreviation ),
				 thickness_scaling = 1.1, right_margin = 11Plots.mm)  #clim =(val_min[k] * 100, val_max[k] * 100)
				Plots.savefig(string("figures\\", folder_name, "\\", metric_name, "_", type))
			end
		end



	elseif model == 3
		folder_name = string(typeof(Liion_electro_chimique()))
		df_ = filter(:couple => ==((E=false,R=false)), df3_result)
		df_E = filter(:couple => ==((E=true,R=false)), df3_result)
		df_R = filter(:couple => ==((E=false,R=true)), df3_result)
		df_ER = filter(:couple => ==((E=true,R=true)), df3_result)

		diff_lifetime = DataFrame(E = (df_.lifetime .- df_E.lifetime)./df_.lifetime, R = (df_.lifetime .- df_R.lifetime)./df_.lifetime, ER= (df_.lifetime .- df_ER.lifetime)./df_.lifetime)
		diff_RES = DataFrame(E = (df_.RES .- df_E.RES)./df_.RES, R = (df_.RES .- df_R.RES)./df_.RES, ER= (df_.RES .- df_ER.RES)./df_.RES)
		diff_cost = DataFrame(E = (df_.cost .- df_E.cost)./df_.cost, R = (df_.cost .- df_R.cost)./df_.cost, ER= (df_.cost .- df_ER.cost)./df_.cost)



		x = [string(unique(df_.pv_size)[i]) for i = 1:length(unique(df_.pv_size))]
		y = [string(unique(df_.bat_size)[i]) for i = 1:length(unique(df_.pv_size))]


		if !isdir(string("figures\\", folder_name))
			mkdir(string("figures\\", folder_name))
		end







		i = 0
		for (metric, data) in enumerate([diff_cost, diff_RES, diff_lifetime])
			i += 1

			if i == 1
				abreviation = "L"
				metric_name = "Lifetime"
			elseif i == 2
				abreviation = "RES"
				metric_name = "RES"
			elseif i == 3
				abreviation = "C"
				metric_name = "Cost"

			end

			for (k,type) in enumerate([:E,:R,:ER])
				z=reshape(data[!,type], length(unique(df_.bat_size)), length(unique(df_.bat_size)))
				Plots.heatmap(x, y, abs.(z*100) , aspect_ratio = 0.7,
				 title = string(metric_name ," impact of ",type," couplage"), xlabel = "PV size", ylabel = "Battery size",
				 colorbar_title = string("\n\n\n(",abreviation, " - ", abreviation, "')/", abreviation ),
				 thickness_scaling = 1.1, right_margin = 11Plots.mm )#clim =(val_min[k] * 100, val_max[k] * 100)
				Plots.savefig(string("figures\\", folder_name, "\\", metric_name, "_", type))
			end
		end



	elseif model == 1
		folder_name = string(typeof(Liion_energy_exchanged()))
		df_ = filter(:couple => ==((E=false,R=false)), df1_result)
		df_E = filter(:couple => ==((E=true,R=false)), df1_result)
		df_R = filter(:couple => ==((E=false,R=true)), df1_result)
		df_ER = filter(:couple => ==((E=true,R=true)), df1_result)

		diff_lifetime = DataFrame(E = (df_.lifetime .- df_E.lifetime)./df_.lifetime, R = (df_.lifetime .- df_R.lifetime)./df_.lifetime, ER= (df_.lifetime .- df_ER.lifetime)./df_.lifetime)
		diff_RES = DataFrame(E = (df_.RES .- df_E.RES)./df_.RES, R = (df_.RES .- df_R.RES)./df_.RES, ER= (df_.RES .- df_ER.RES)./df_.RES)
		diff_cost = DataFrame(E = (df_.cost .- df_E.cost)./df_.cost, R = (df_.cost .- df_R.cost)./df_.cost, ER= (df_.cost .- df_ER.cost)./df_.cost)



		x = [string(unique(df_.pv_size)[i]) for i = 1:length(unique(df_.pv_size))]
		y = [string(unique(df_.bat_size)[i]) for i = 1:length(unique(df_.pv_size))]


		if !isdir(string("figures\\", folder_name))
			mkdir(string("figures\\", folder_name))
		end







		i = 0
		for (metric, data) in enumerate([diff_cost, diff_RES, diff_lifetime])
			i += 1

			if i == 1
				abreviation = "L"
				metric_name = "Lifetime"
			elseif i == 2
				abreviation = "RES"
				metric_name = "RES"
			elseif i == 3
				abreviation = "C"
				metric_name = "Cost"

			end

			for (k,type) in enumerate([:E,:R,:ER])
				z=reshape(data[!,type], length(unique(df_.bat_size)), length(unique(df_.bat_size)))
				Plots.heatmap(x, y, abs.(z*100) , aspect_ratio = 0.7,
				 title = string(metric_name ," impact of ",type," couplage"), xlabel = "PV size", ylabel = "Battery size",
				 colorbar_title = string("\n\n\n(",abreviation, " - ", abreviation, "')/", abreviation ),
				 thickness_scaling = 1.1, right_margin = 11Plots.mm) #clim =(val_min[k] * 100, val_max[k] * 100)
				Plots.savefig(string("figures\\", folder_name, "\\", metric_name, "_", type))
			end
		end


	end

end










using StatsPlots



anim = @animate for i in 5:5:150


	df = filter(:pv_size => ==(i), df_result)


	@df df Plots.scatter(
	    :bat_size,
	    :lifetime,
	    group = :couple,
	    title = string("pv_size = ",i),
	    xlabel = "bat_size",
	    ylabel = "life time",
	    m = (0.5, [:circle :vline :hline :+], 12),
	    bg = RGB(0.2, 0.2, 0.2)
	)
end

gif(anim, "anim_fps.gif", fps = 1)
