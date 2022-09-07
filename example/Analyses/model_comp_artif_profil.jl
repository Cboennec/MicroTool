
# Load packages
include("..\\..\\src\\Genesys.jl")

using Main.Genesys
using JLD, Dates, Seaborn
using CSV, DataFrames


pygui(true)

# Parameters of the simulation
const nh, ny, ns = 8760, 50, 1

# Load input data
data = load(joinpath("example","data","ausgrid_5_twostage.jld"))

#https://www.researchgate.net/figure/Battery-cycle-to-failure-versus-the-DOD-gap-34_fig5_276136726
fatigue_data = DataFrame(CSV.File("example\\data\\fatigue_data_NMC.csv", delim = ",", header = [Symbol("DoD"),Symbol("cycle")], types=Dict(:DoD=>Float64, :cycle=>Float64)))



#####################################################################""
#Test part





#############################################

soh_treshold = 0.8

# 1/j et 2/J versions pour 12 outputs
# maxi = [1., .5, 1., 1., .2, .8, 1., .5, 1., 1., .2, .8]
# mini = [.0, 0., .5, .8, .0, .2, .0, .0, .5, .8, .0, .2]
# half_wave_length = [repeat([12],6); repeat([6],6),]

nb_freq = 2
nb_ampl = 6

maxi = repeat([1., 1., .5, 1., .2, .8],nb_freq)
mini = repeat([.0, .5, 0., .8, .0, .2],nb_freq)

half_wave_length = [repeat([12],nb_ampl); repeat([6],nb_ampl)]#; repeat([1],nb_ampl)]


nb_profil = length(maxi)

lin_res = []
ec_res = []
rf_res = []



for i in 1:nb_profil

	croiss = [(h)/half_wave_length[i] for h in 1:(half_wave_length[i]-1)]
	rise = [0; croiss]
	decroiss = reverse(croiss)
	fall = [1; decroiss]

	oneday = repeat([rise; fall].* (maxi[i]-mini[i]) .+ mini[i], Int(12/half_wave_length[i]))

	artif = repeat(oneday, 365, ny+1)
	artif = vcat(artif, transpose(artif[end,:]))







	microgrid_rainflow = Microgrid(parameters = GlobalParameters(nh, ny, ns, renewable_share = 1.))
	add!(microgrid_rainflow, Demand(carrier = Electricity()),
				Solar(),
				Liion_rainflow(soc_model = "artificial",couplage = (E=false,R=false), artificial_soc_profil = artif, fatigue_data = fatigue_data, calendar = true),
				Grid(carrier = Electricity()))




	ω_e = Scenarios(microgrid_rainflow, data["ω_optim"], adjust_length=true)


		controller = RBC(options = RBCOptions(policy_selection =  2))
		designer = Manual(generations = [10.], storages = [10.], subscribed_power = [5.])

		designer_ec = initialize_designer!(microgrid_rainflow, designer, ω_e)
		controller_ec = initialize_controller!(microgrid_rainflow, controller, ω_e)

		simulate!(microgrid_rainflow, controller_ec, designer_ec, ω_e, options = Genesys.Options(mode = "serial"))





	microgrid_ec = Microgrid(parameters = GlobalParameters(nh, ny, ns, renewable_share = 1.))
	add!(microgrid_ec, Demand(carrier = Electricity()),
	                Solar(),
	                Liion_electro_chimique(soc_model = "artificial", couplage = (E=false,R=false), artificial_soc_profil = artif),
	                Grid(carrier = Electricity()))

	ω_e = Scenarios(microgrid_ec, data["ω_optim"], adjust_length=true)


	controller = RBC(options = RBCOptions(policy_selection =  2))
	designer = Manual(generations = [10.], storages = [10.], subscribed_power = [5.])

	designer_ec = initialize_designer!(microgrid_ec, designer, ω_e)
	controller_ec = initialize_controller!(microgrid_ec, controller, ω_e)

	simulate!(microgrid_ec, controller_ec, designer_ec, ω_e, options = Genesys.Options(mode = "serial"))





	microgrid_lin = Microgrid(parameters = GlobalParameters(nh, ny, ns, renewable_share = 1.))
	add!(microgrid_lin, Demand(carrier = Electricity()),
	                Solar(),
	                Liion_energy_exchanged(soc_model = "artificial", couplage = (E=false,R=false), artificial_soc_profil = artif, calendar = true, nCycle = fatigue_data.cycle[findfirst(fatigue_data.DoD .> (0.6))]),
	                Grid(carrier = Electricity()))

	ω_e = Scenarios(microgrid_lin, data["ω_optim"], adjust_length=true)


	controller = RBC(options = RBCOptions(policy_selection =  2))
	designer = Manual(generations = [10.], storages = [10.], subscribed_power = [5.])

	designer_ec = initialize_designer!(microgrid_lin, designer, ω_e)
	controller_ec = initialize_controller!(microgrid_lin, controller, ω_e)

	simulate!(microgrid_lin, controller_ec, designer_ec, ω_e, options = Genesys.Options(mode = "serial"))




	id = findfirst((microgrid_lin.storages[1].soh[:,2:(end-1),1]) .< soh_treshold)
	#if nothing, then it mean that the threshold has not been reached yet
	if isnothing(id)
		#then we extrapolate from the begining for linear models
		push!(lin_res, (1- soh_treshold) / (1-(microgrid_lin.storages[1].soh[end-1,(end-1),1])) * ny )
		println("lin did not finish for profil : ",i)
	else
		push!(lin_res, id[2] + id[1]/8760)
	end


	id = findfirst((microgrid_ec.storages[1].soh[:,2:(end-1),1]) .< soh_treshold)
	if isnothing(id)
		#And we extrapolate from the last year for non linear models since ec tend to be linear after some years
		(1- soh_treshold) / (1-(microgrid_ec.storages[1].soh[end-1,(end-1),1])) * ny
		last_year = (1-(microgrid_ec.storages[1].soh[end-1,(end-2),1])) - (1-(microgrid_ec.storages[1].soh[end-1,(end-1),1]))
		rest = ((1-(microgrid_ec.storages[1].soh[end-1,(end-1),1])) - soh_treshold) / last_year
		push!(ec_res, ny + rest)
		println("ec did not finish for profil : ",i)

	else
		push!(ec_res, id[2] + id[1]/8760)
	end

	id = findfirst((microgrid_rainflow.storages[1].soh[:,2:(end-1),1]) .< soh_treshold)
	if isnothing(id)
		#In this case rainflow is considered linear since the cycling is fully homogenous
		push!(rf_res, (1- soh_treshold) / (1-(microgrid_rainflow.storages[1].soh[1,(end),1])) * ny)
		println("rf did not finish for profil : ",i)
	else
		push!(rf_res, id[2] + id[1]/8760)
	end

end


using Pandas


profil_num = repeat(1:nb_profil, 3)
model_name = [repeat(["Exchanged Energy"], nb_profil); repeat(["Rainflow"], nb_profil); repeat(["Electro-Chemical"], nb_profil) ]
delta_SOH = [lin_res ; rf_res; ec_res]


df_pandas = Pandas.DataFrame(Dict("delta_SOH" => delta_SOH, "profile_num" => profil_num, "model_name" => model_name) )


figure("Model comparison on artificial profiles")
Seaborn.set_context("paper", rc=Dict("axes.labelsize"=>28, "legend.fontsize"=> 22,  "legend.title_fontsize"=> 28,
"ytick.labelsize" => 22, "xtick.labelsize" => 22))

ax  = Seaborn.catplot(x="profile_num", y="delta_SOH",
			   hue="model_name",
			   data=df_pandas, kind="bar",
			   palette="Set2",
			   height=4, aspect=2.5, legend = false);

Seaborn.legend( loc="upper left", title ="Ageing Models")

ax.set_ylabels("Time before EoL (years)")
ax.set_xlabels("SoC profiles")

for i in 1:(nb_freq-1)
	axvline(x=(nb_ampl * i) - 0.5, ymin=0, ymax=1, color="black")
end
