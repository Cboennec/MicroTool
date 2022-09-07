
# Load packages
include("..\\src\\Genesys.jl")

using Main.Genesys
using JLD, Dates, Seaborn
using CSV, DataFrames


pygui(true)

# Parameters of the simulation
const nh, ny, ns = 8760, 15, 1

# Load input data
data = load(joinpath("example","data","ausgrid_5_twostage.jld"))

#https://www.researchgate.net/figure/Battery-cycle-to-failure-versus-the-DOD-gap-34_fig5_276136726
fatigue_data = DataFrames.DataFrame(CSV.File("example\\data\\fatigue_data_NMC.csv", delim = ",", header = [Symbol("DoD"),Symbol("cycle")], types=Dict(:DoD=>Float64, :cycle=>Float64)))











mg = Microgrid(parameters = GlobalParameters(nh, ny, ns, renewable_share = 1.))

add!(mg, Demand(carrier = Electricity()),
		  Solar(),
		  Liion_electro_chimique(update_by_year = 12, soc_model = "linear", couplage= (E=true, R=true)),
		  Grid(carrier = Electricity()))


ω_e = Scenarios(mg, data["ω_optim"], adjust_length=true)


	controller = RBC(options = RBCOptions(policy_selection =  2))
	designer = Manual(generations = [10.], storages = [10.], subscribed_power = [5.])

	designer_ec = initialize_designer!(mg, designer, ω_e)
	controller_ec = initialize_controller!(mg, controller, ω_e)



	simulate!(mg, controller_ec, designer_ec, ω_e, options = Genesys.Options(mode = "serial"))

	metrics = Metrics(mg, designer_ec)

	# Plots
	#plot_operation(mg, y=2:ny)



efficiency = (mg.storages[1].soc[2:end-1,1:end-1,1] .- mg.storages[1].soc[3:end,1:end-1,1]) .* mg.storages[1].soh[3:end,1:end-1,1] ./ controller_ec.decisions.storages[1][2:end,:,1]

discharging_index = (mg.storages[1].soc[2:end-1,1:end-1,1] .- mg.storages[1].soc[3:end,1:end-1,1]) .> 0

(mg.storages[1].soc[2:end-1,1:end-1,1] .- mg.storages[1].soc[3:end,1:end-1,1])
charging_index1 = (mg.storages[1].soc[2:end-1,1:end-1,1] .- mg.storages[1].soc[3:end,1:end-1,1]) .< 0
charging_index2 =  controller_ec.decisions.storages[1][2:end,:,1] .< 0
charging_index = charging_index1 .& charging_index2

plot(efficiency[charging_index])

plot(1 ./ efficiency[discharging_index])





function compute_operation_soc_tremblay_dessaint(liion::AbstractLiion, state::NamedTuple{(:Erated, :soc, :soh), Tuple{Float64, Float64, Float64}}, V::Float64, decision::Float64, Δh::Int64)

	# Control power constraint and correction
	E0 = liion.tremblay_dessaint_params.E0
	A = liion.tremblay_dessaint_params.A
	B = liion.tremblay_dessaint_params.B
	K = liion.tremblay_dessaint_params.K

	# up to 70% boost on resistance linearly increase from ini to battery end of life here at 70% soh
	if liion.couplage.R
	   R = liion.tremblay_dessaint_params.R * (1 + (0.7 * (1-state.soh) / (1-(liion.SoH_threshold))))
	else
	   R = liion.tremblay_dessaint_params.R
	end

	I_max = liion.tremblay_dessaint_params.I_max
	soc_min = liion.α_soc_min
	soc_max = liion.α_soc_max

	Qnom = state.Erated / (V * liion.tremblay_dessaint_params.Npara * liion.tremblay_dessaint_params.Nserie) # TODO 3.7 en V pour convertir Erated de kWh en Ah
	if liion.couplage.E
	   Qnom *= state.soh #TODO faire verifier
	end


	I = decision / (V * liion.tremblay_dessaint_params.Npara * liion.tremblay_dessaint_params.Nserie)

	it = (1-state.soc)*Qnom

	# calcul du soc suivant et volt suivant

	it_next = it + I * Δh
	soc_next = 1 - it_next/Qnom



	soc_next = min(max(soc_next,soc_min),soc_max)
	#it_next = (1-soc_next)*Qnom
	#I = (it_next - it)/Δh


	#decision = I * V
	power_dch = max(min(decision, V*I_max ), 0.)
	power_ch = min(max(decision,-V*I_max ), 0.)
	I = min(max(-I_max,I),I_max)


	if I>0
		# volt_next = E0 - R*I - K*Qnom * (1 - soc_next) / soc_next - K*I / soc_next + A*exp(-B * Qnom * (1 - soc_next))
		volt_next = E0 -R*I - K*it_next * Qnom /(Qnom - it_next) - K*I * Qnom /(Qnom - it_next) + A*exp(-B * it_next)
	else
		# volt_next = E0 -R*I - K*Qnom * (1 - soc_next) / soc_next - K*I / (0.9 - soc_next) + A*exp(-B * Qnom * (1 - soc_next))
		volt_next = E0 -R*I - K*it_next * Qnom /(Qnom - it_next) - K*I * Qnom /(it_next - 0.1*Qnom) + A*exp.(-B * it_next)
	end
	#volt_next > 4 ? volt_next = 4 : nothing


	return soc_next, volt_next, (power_ch + power_dch)* liion.tremblay_dessaint_params.Npara* liion.tremblay_dessaint_params.Nserie, I
end
