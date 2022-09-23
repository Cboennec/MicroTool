#=
    Li-ion battery modelling
 =#

#ref :  Optimal Sizing and Control of a PV-EV-BES Charging System Including Primary Frequency Control and Component Degradation
# Wiljan Vermeer

#original paper : Degradation of lithium ion batteries employing graphite negatives and nickelecobaltemanganese oxide þ spinel manganese oxide positives: Part 1, aging mechanisms and life estimation
#John Wang


#fit for c1
# using Polynomials
# # Here the values of the paper  https://doi.org/10.1016/j.jpowsour.2014.07.030 using webplot digitilizer
# points = [283.1958762886598 0.002115384615384617;
# 293.09278350515467 0.00078846153846154;
# 307.0103092783505 0.0010000000000000009;
# 319.07216494845363 0.0045000000000000005]
# x = points[:,1]
# y = points[:,2]
# quadfit=Polynomials.fit(x,y,2)
#
# using Plots
#
# Plots.plot(x,y,label="Data")
#
# plot!(quadfit,x[1],x[end],label="Quadratic Fit")


soc_model_names = ["tremblay_dessaint", "linear", "vermeer", "artificial"]

mutable struct Vermeer_params
   c1::Float64
   c2::Float64
   c3::Float64
   c4::Float64
   c6::Float64
   Vermeer_params(; c1 = 0.0004111236450494715 , # 8.699477998200086e-6 * T² - 0.005177224156155348 * T + 0.7706741325793118 ; T=298.15
    c2 = 3.524e-1, #-6.7e-3 * T + 2.35
    c3 = 50,
    c4 = 1.035,
    c6 = 0.0000022831 #2% / 8760
	) = new( c1, c2, c3, c4, c6)
end




   mutable struct Liion_vermeer <: AbstractLiion

	# Parameters
    α_p_ch::Float64
   	α_p_dch::Float64
    η_ch::Float64 #Charging yield / efficacity
   	η_dch::Float64 #Discharging yield / efficacity
   	η_self::Float64 #Auto discarge factor
   	α_soc_min::Float64 #min threshold of charge (normalized)
   	α_soc_max::Float64 #max threshold of charge (normalized)
   	lifetime::Int64
   	nCycle::Float64
   	bounds::NamedTuple{(:lb, :ub), Tuple{Float64, Float64}}
   	SoH_threshold::Float64 # SoH level to replace battery
   	couplage::NamedTuple{(:E, :R), Tuple{Bool, Bool}}  #a boolean tuple to tell wether or not the soh should influence the other parameters.
	Npara::Int64
	Nseries::Int64

	#needed for SoH computation but currently used as a constant
	temperature::Float64

   	#Model dynamics
   	soc_model::String #model name

   	# Initial conditions
   	Erated_ini::Float64  # capacité de la batterie en Wh
   	soc_ini::Float64 # first state of charge for the begining of simulation
   	soh_ini::Float64 # first state of health for the begining of simulation

   	# Variables
	ΔE_tot::AbstractArray{Float64,1} # This is used as a memory of the cumulated fatigue of the battery, see ref [1] for some détails
   	Erated::AbstractArray{Float64,2} # Battery capacity
   	carrier::EnergyCarrier #Type of energy
   	soc::AbstractArray{Float64,3} #3 dim matrix (h,y,s) containing the state of charge [0-1]
   	soh::AbstractArray{Float64,3} #3 dim matrix (h,y,s) containing the state of health [0-1]
	voltage ::Array{Float64,3} # evolution de la tension

   	#tremblay_dessaint
   	tremblay_dessaint_params::Tremblay_dessaint_params # A list of params only use for tremblay dessaint
	vermeer_params::Vermeer_params
   	# Eco
   	cost::AbstractArray{Float64,2}

   	# Inner constructor
   	Liion_vermeer(; α_p_ch = 1.5,
   		α_p_dch = 1.5,
   		η_ch = 0.9,
   		η_dch = 0.9,
   		η_self = 0.0005,
   		α_soc_min = 0.2,
   		α_soc_max = 0.8,
   		lifetime = 12,
   		nCycle = 2500.,
   		bounds = (lb = 0., ub = 1000.),
   		SoH_threshold = 0.8,
   		couplage = (E = false, R = false),
		Npara = 1,
		Nseries = 1,
		temperature = 293,
   		soc_model = "linear",
   		Erated_ini = 1e-6,
   		soc_ini = 0.5,
   		soh_ini = 1.,
   		update_by_year = 1) =  verification_liion_params(α_p_ch, α_p_dch, η_ch, η_dch, η_self, α_soc_min, α_soc_max, lifetime, nCycle, bounds,
   			SoH_threshold, couplage, Npara, Nseries, soc_model, Erated_ini, soc_ini, soh_ini) ?
   			new(α_p_ch, α_p_dch, η_ch, η_dch, η_self, α_soc_min, α_soc_max, lifetime, nCycle, bounds,
   			SoH_threshold, couplage, Npara, Nseries, temperature, soc_model, Erated_ini, soc_ini, soh_ini) : nothing

end

### Preallocation
 function preallocate!(liion::Liion_vermeer, nh::Int64, ny::Int64, ns::Int64)
     liion.Erated = convert(SharedArray,zeros(ny+1, ns)) ; liion.Erated[1,:] .= liion.Erated_ini
     liion.carrier = Electricity()
     liion.carrier.power = convert(SharedArray,zeros(nh, ny, ns))
     liion.soc = convert(SharedArray,zeros(nh+1, ny+1, ns)) ; liion.soc[1,1,:] .= liion.soc_ini
     liion.soh = convert(SharedArray,zeros(nh+1, ny+1, ns)) ; liion.soh[1,1,:] .= liion.soh_ini
     liion.cost = convert(SharedArray,zeros(ny, ns))
	 liion.tremblay_dessaint_params = Tremblay_dessaint_params()
	 liion.voltage = convert(SharedArray,zeros(nh+1, ny+1, ns)) ; liion.voltage[1,1:2,:] .= 3.7
	 liion.vermeer_params = Vermeer_params()
	 liion.ΔE_tot = convert(SharedArray,zeros(ns))

     return liion
 end

 ### Operation dynamic
function compute_operation_dynamics!(h::Int64, y::Int64, s::Int64, liion::Liion_vermeer, decision::Float64, Δh::Int64)

	if liion.soc_model == "vermeer"
		liion.soc[h+1,y,s], liion.carrier.power[h,y,s] = compute_operation_soc_Vermeer(liion, (Erated = liion.Erated[y,s], soc = liion.soc[h,y,s], soh = liion.soh[h,y,s]), decision, Δh)
	elseif liion.soc_model == "tremblay_dessaint"
		liion.soc[h+1,y,s], liion.voltage[h+1,y,s], liion.carrier.power[h,y,s], liion.current[h,y,s] = compute_operation_soc_tremblay_dessaint(liion, (Erated = liion.Erated[y,s], soc = liion.soc[h,y,s], soh = liion.soh[h,y,s]),  liion.voltage[h,y,s], decision, Δh)
	elseif  liion.soc_model == "linear"
		liion.soc[h+1,y,s], liion.carrier.power[h,y,s]  = compute_operation_soc_linear(liion, (Erated = liion.Erated[y,s], soc = liion.soc[h,y,s], soh = liion.soh[h,y,s]), decision, Δh)
	end

	liion.soh[h+1,y,s] = liion.soh[h,y,s] - (0.1/8760)
	#liion.soh[h+1,y,s], liion.ΔE_tot[s] = compute_operation_soh_vermeer(liion, (Erated = liion.Erated[y,s], soc = liion.soc[h,y,s], soh = liion.soh[h,y,s]), decision, Δh, liion.soc[h+1,y,s], liion.ΔE_tot[s])

end




 ########################
 ##### SOH models #######
 ########################


 function compute_operation_soh_vermeer(liion::Liion_vermeer, state::NamedTuple{(:Erated, :soc, :soh), Tuple{Float64, Float64, Float64}}, decision::Float64, Δh::Int64, soc_next::Float64, ΔE_tot::Float64 )

	 η_ini = 0.98   #Fixed (dis)charging efficiency for both BES and EV (0.98)     dans la nomenclature
	 η = η_ini - ((1-state.soh)/12)   #(15) simplifié

	power_dch = max(min(decision, liion.α_p_dch * state.Erated, state.soh *  state.Erated / Δh, η * (liion.α_soc_max - state.soc) * state.Erated * state.soh / Δh), 0.)
  	power_ch = min(max(decision, -liion.α_p_ch * state.Erated, -state.soh *  state.Erated / Δh, (liion.α_soc_min - state.soc) * state.Erated * state.soh / Δh / η), 0.)

	Voc_linear = liion.Nseries * (3.42 + 0.7 * state.soc)

    P = (η * power_ch) + (power_dch / η) #(13)

	I_cell =  P/(liion.Npara * Voc_linear)

	DoD = abs(soc_next - state.soc)

	ΔE_percent = liion.vermeer_params.c1 * exp(liion.vermeer_params.c2 * abs(I_cell)) * ((liion.vermeer_params.c3 * DoD) / liion.vermeer_params.c4) * abs(I_cell) * Δh

	ΔE_cycle =  ΔE_percent * state.Erated / 100

	ΔE_cal = liion.vermeer_params.c6 * Δh * state.Erated

	#println("cal = ", ΔE_cal , " cycle = ", ΔE_cycle)
	ΔE_tot += (ΔE_cycle + ΔE_cal)

	return 1 - (ΔE_tot /state.Erated) , ΔE_tot
end

 ### Investment dynamic
 function compute_investment_dynamics!(y::Int64, s::Int64, liion::Liion_vermeer, decision::Union{Float64, Int64})
	 liion.Erated[y+1,s], liion.soc[1,y+1,s], liion.soh[1,y+1,s], liion.voltage[1,y+1,s], liion.ΔE_tot[s] = compute_investment_dynamics(liion, (Erated = liion.Erated[y,s], soc = liion.soc[end,y,s], soh = liion.soh[end,y,s], voltage = liion.voltage[end,y,s], ΔE_tot = liion.ΔE_tot[s]), decision)
 end

 function compute_investment_dynamics(liion::Liion_vermeer, state::NamedTuple{(:Erated, :soc, :soh, :voltage, :ΔE_tot), Tuple{Float64, Float64, Float64, Float64, Float64}}, decision::Union{Float64, Int64})
	 if decision > 1e-2
		 Erated_next = decision
		 soc_next = liion.soc_ini
		 soh_next =  1.
		 ΔE_tot_next = 0.

	 else
		 Erated_next = state.Erated
		 soc_next = state.soc
		 soh_next = state.soh
		 ΔE_tot_next = state.ΔE_tot

	 end

	#TODO faire verifier
	if state.voltage != 0
		voltage_next = state.voltage
	else
	   voltage_next = liion.soc_ini*0.7 + 3.3
	end

	 return Erated_next, soc_next, soh_next, voltage_next, ΔE_tot_next
 end


#for rule based #TODO change computataion for soc and soh
 function compute_operation_dynamics(liion::Liion_vermeer, state::NamedTuple{(:Erated, :soc, :soh), Tuple{Float64, Float64, Float64}}, decision::Float64, Δh::Int64)
      # Control power constraint and correction
	  if liion.couplage.E
     	 Erated = state.Erated * state.soh
      else
     	 Erated = state.Erated
      end

	  power_dch = max(min(decision, liion.α_p_dch * Erated, state.soh * Erated / Δh, liion.η_dch * (state.soc * (1. - liion.η_self * Δh) - liion.α_soc_min) * Erated / Δh), 0.)
      power_ch = min(max(decision, -liion.α_p_ch * Erated, -state.soh * Erated / Δh, (state.soc * (1. - liion.η_self * Δh) - liion.α_soc_max) * Erated / Δh / liion.η_ch), 0.)

	  #for rule base, doesnt matter
	  soc_next = state.soc * (1. - liion.η_self * Δh) - (power_ch * liion.η_ch + power_dch / liion.η_dch) * Δh / state.Erated
	  # SoH dynamic
	  soh_next = state.soh - (power_dch - power_ch) * Δh / (2. * liion.nCycle * (liion.α_soc_max - liion.α_soc_min) * state.Erated)
      return soc_next, soh_next, power_dch + power_ch
 end



 function verification_liion_params(α_p_ch::Float64, α_p_dch::Float64, η_ch::Float64, η_dch::Float64, η_self::Float64,
 	α_soc_min::Float64, α_soc_max::Float64, lifetime::Int64, nCycle::Float64, bounds::NamedTuple{(:lb, :ub), Tuple{Float64, Float64}},
 	SoH_threshold::Float64, couplage::NamedTuple{(:E,:R), Tuple{Bool,Bool}}, Npara::Int64, Nseries::Int64, soc_model::String, Erated_ini::Float64, soc_ini::Float64,
 	soh_ini::Float64)

 	validation = true

 	# α_p_ch, α_p_dch
 	if bounds.lb > bounds.ub
 		error("Lower bound cannot be inferieur to upper bound")
 		validation = false
 	end

 	if α_soc_min > α_soc_max
 		error("Lower bound cannot be inferieur to upper bound (soc bounds)")
 		validation = false
 	end

 	if !(soc_model in soc_model_names)
 		error(soc_model ," is not an authorized Liion state of charge model. you need to pick one from the following list : ", soc_model_names)
 		validation = false
 	end

 	if η_ch < 0 || η_ch > 1
 		error("η_ch out of bound [0-1] with value :" , η_ch)
 		validation = false
 	end

 	if η_dch < 0 || η_dch > 1
 		error("η_dch out of bound [0-1] with value :" , η_dch)
 		validation = false
 	end

 	if η_self < 0 || η_self > 1
 		error("η_self out of bound [0-1] with value :" , η_self)
 		validation = false
 	end

 	if α_soc_min < 0 || α_soc_min > 1
 		error("α_soc_min out of bound [0-1] with value :" , α_soc_min)
 		validation = false
 	end

 	if α_soc_max < 0 || α_soc_max > 1
 		error("α_soc_max out of bound [0-1] with value :" , α_soc_max)
 		validation = false
 	end

 	if Erated_ini < 0
 		error("Erated_ini must be positive")
 		validation = false
 	end

 	if soc_ini < 0 || soc_ini > 1
 		error("soc_ini out of bound [0-1] with value :" , soc_ini)
 		validation = false
 	end

 	if soh_ini < 0 || soh_ini > 1
 		error("soh_ini out of bound [0-1] with value :" , soh_ini)
 		validation = false
 	end

 	if SoH_threshold  < 0 || SoH_threshold >= 1
 		error("SoH_threshold out of bound [0-1[ with value :" , SoH_threshold)
 		validation = false
 	end

 	return validation
 end
