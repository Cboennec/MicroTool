#=
    Li-ion battery modelling
 =#

#[1] Modeling of Lithium-Ion Battery Degradation for Cell Life Assessment
#Authors : Bolun Xu, Student Member, IEEE, Alexandre Oudalov, Andreas Ulbig, Member, IEEE, Göran Andersson, Fellow, IEEE, and Daniel S. Kirschen, Fellow, IEEE

soc_model_names = ["tremblay_dessaint", "linear", "vermeer", "artificial"]

mutable struct Electro_chimique_params
   alpha_sei::Float64
   beta_sei::Float64
   k_delta1::Float64
   k_delta2::Float64
   k_sigma::Float64
   sigma_ref::Float64
   k_T::Float64
   T_ref::Float64
   k_t::Float64



   # for NMC parameters #
   #Technoeconomic model of second-life batteries for utility-scale solar
   #considering calendar and cycle aging
   #Ian Mathewsa,⁎, Bolun Xub, Wei Hea, Vanessa Barretoc, Tonio Buonassisia, Ian Marius Peters

   Electro_chimique_params(; alpha_sei = 5.75e-2,
    beta_sei = 121,
    k_delta1 = 1.0487e-4,
    k_delta2 = 2.03,
    k_sigma = 1.04,
    sigma_ref = 0.5,
   	k_T = 6.93e-2,
   	T_ref = 298,
   	k_t = 4.14e-10) = new( alpha_sei, beta_sei, k_delta1, k_delta2, k_sigma, sigma_ref, k_T, T_ref, k_t)
end



#rainflow ref : Optimal Battery Control Under Cycle Aging Mechanisms in Pay for Performance Settings
# Yuanyuan Shi, Bolun Xu, Yushi Tan, Daniel Kirschen, Baosen Zhang
   mutable struct Liion_electro_chimique <: AbstractLiion

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

	#needed for SoH computation but currently used as a constant
	temperature::Float64

   	#Model dynamics
   	soc_model::String #model name

   	# Initial conditions
   	Erated_ini::Float64  # capacité de la batterie en Wh
   	soc_ini::Float64 # first state of charge for the begining of simulation
   	soh_ini::Float64 # first state of health for the begining of simulation

   	#input specific params
   	update_by_year::Int64 # Rainflow Soh computation by year

	artificial_soc_profil::Array{Float64,2} # soc profil for study purpose

   	# Variables
	Sum_fd::AbstractArray{Float64,1} # This is used as a memory of the cumulated fatigue of the battery, see ref [1] for some détails
   	Erated::AbstractArray{Float64,2} # Battery capacity
   	carrier::EnergyCarrier #Type of energy
   	soc::AbstractArray{Float64,3} #3 dim matrix (h,y,s) containing the state of charge [0-1]
   	soh::AbstractArray{Float64,3} #3 dim matrix (h,y,s) containing the state of health [0-1]
   	#tremblay_dessaint
   	voltage ::Array{Float64,3} # evolution de la tension
   	current ::Array{Float64,3} # evolution du courant ?
   	tremblay_dessaint_params::Tremblay_dessaint_params # A list of params only use for tremblay dessaint
	Electro_chimique_params::Electro_chimique_params # A list of spécifique params for this type of SoH computation

   	# Eco
   	cost::AbstractArray{Float64,2}

   	# Inner constructor
   	Liion_electro_chimique(; α_p_ch = 1.5,
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
   		couplage = (E = true, R = false),
		temperature = 298,
   		soc_model = "linear",
   		Erated_ini = 1e-6,
   		soc_ini = 0.5,
   		soh_ini = 1.,
   		update_by_year = 12,
		artificial_soc_profil = zeros(8760,1)) =  verification_liion_params(α_p_ch, α_p_dch, η_ch, η_dch, η_self, α_soc_min, α_soc_max, lifetime, nCycle, bounds,
   			SoH_threshold, couplage, soc_model, Erated_ini, soc_ini, soh_ini, update_by_year, artificial_soc_profil) ?
   			new(α_p_ch, α_p_dch, η_ch, η_dch, η_self, α_soc_min, α_soc_max, lifetime, nCycle, bounds,
   			SoH_threshold, couplage, temperature, soc_model, Erated_ini, soc_ini, soh_ini, update_by_year, artificial_soc_profil) : nothing

end

### Preallocation
 function preallocate!(liion::Liion_electro_chimique, nh::Int64, ny::Int64, ns::Int64)
     liion.Erated = convert(SharedArray,zeros(ny+1, ns)) ; liion.Erated[1,:] .= liion.Erated_ini
     liion.carrier = Electricity()
     liion.carrier.power = convert(SharedArray,zeros(nh, ny, ns))
	 if liion.soc_model == "artificial"
		 liion.soc = convert(SharedArray,reshape(repeat(liion.artificial_soc_profil,ns), (nh+1,ny+1,ns)))
	 else
		 liion.soc = convert(SharedArray,zeros(nh+1, ny+1, ns)) ; liion.soc[1,1,:] .= liion.soc_ini
     end
     liion.soh = convert(SharedArray,zeros(nh+1, ny+1, ns)) ; liion.soh[1,1,:] .= liion.soh_ini
	 liion.voltage = convert(SharedArray,zeros(nh+1, ny+1, ns)) ; liion.voltage[1,1:2,:] .= 3.7
	 liion.current = convert(SharedArray,zeros(nh+1, ny+1, ns))
     liion.cost = convert(SharedArray,zeros(ny, ns))
	 liion.tremblay_dessaint_params = Tremblay_dessaint_params()
	 liion.Electro_chimique_params = Electro_chimique_params()
	 liion.Sum_fd = convert(SharedArray,zeros(ns))
     return liion
 end

 ### Operation dynamic
function compute_operation_dynamics!(h::Int64, y::Int64, s::Int64, liion::Liion_electro_chimique, decision::Float64, Δh::Int64)

	if liion.soc_model == "tremblay_dessaint"
		liion.soc[h+1,y,s], liion.voltage[h+1,y,s], liion.carrier.power[h,y,s], liion.current[h,y,s] = compute_operation_soc_tremblay_dessaint(liion, (Erated = liion.Erated[y,s], soc = liion.soc[h,y,s], soh = liion.soh[h,y,s]),  liion.voltage[h,y,s], decision, Δh)
	elseif  liion.soc_model == "linear"
		liion.soc[h+1,y,s], liion.carrier.power[h,y,s] = compute_operation_soc_linear(liion, (Erated = liion.Erated[y,s], soc = liion.soc[h,y,s], soh = liion.soh[h,y,s]), decision, Δh)
	elseif liion.soc_model == "artificial"

	end


	h_between_update = convert(Int64,floor(8760/liion.update_by_year))
	#SoH computation
    if (h%h_between_update) != 0
		liion.soh[h+1,y,s] = liion.soh[h,y,s]
    else #rainflow computaion
		interval = (h-h_between_update+1):h

		if isnan(liion.soc[interval[2],y,s])
			println("y = : ",y, ";     soc = ", liion.soc[h,y,s])
			println("interval = : ", interval)
		end

		liion.soh[h+1,y,s], liion.Sum_fd[s] = compute_operation_soh_rainflow(liion, (Erated = liion.Erated[y,s], soc = liion.soc[h,y,s], soh = liion.soh[h,y,s]), decision, Δh,  liion.soc[interval,y,s], liion.Sum_fd[s])
    end
end




 ########################
 ##### SOH models #######
 ########################


 function compute_operation_soh_rainflow(liion::Liion_electro_chimique, state::NamedTuple{(:Erated, :soc, :soh), Tuple{Float64, Float64, Float64}}, decision::Float64, Δh::Int64, soc::Vector{Float64}, Sum_fd::Float64)

	soc_peak, soc_peak_id = get_soc_peaks(soc)
	#Then compute the DoD sequence by extracting the subcycles DoD

	DoD_seq = Float64[] #Sequence of all the charging and decharging half cycles DoDs
	mean_Soc_seq = Float64[]
	delta_t_seq = Int64[]

	i = 1


	while i+3 <= length(soc_peak_id)	#Define your 3 deltas with 4 consecutives points
		delta1 = abs( soc[soc_peak_id[i+1]] - soc[soc_peak_id[i]] )
	 	delta2 = abs( soc[soc_peak_id[i+2]] - soc[soc_peak_id[i+1]] )
	 	delta3 = abs( soc[soc_peak_id[i+3]] - soc[soc_peak_id[i+2]] )

	 	#rainflow sub-cycle criterion
	 	if delta2 <= delta1 && delta2 <= delta3
			push!(DoD_seq, delta2) #1 half cycle of DoD delta2 +
			push!(DoD_seq, delta2) #1 half cycle
			push!(mean_Soc_seq, (soc[soc_peak_id[i+2]] + soc[soc_peak_id[i+1]]) /2 ) #SoC mean
			push!(mean_Soc_seq, (soc[soc_peak_id[i+2]] + soc[soc_peak_id[i+1]]) /2 ) #SoC mean
			push!(delta_t_seq, soc_peak_id[i+2] - soc_peak_id[i+1])#delta_t
			push!(delta_t_seq, soc_peak_id[i+2] - soc_peak_id[i+1])#delta_t, we use the time of the second half cycle because the first one is altered by the algorithme.


	 		deleteat!(soc_peak_id, i+2) #start with i+2 index or you will delete i+1 and i+3
	 		deleteat!(soc_peak_id, i+1)
		else #else use the following point sequence
	 		i = i+1
	 	end
 	end

 	#Then add the englobing (those who make the other cycles "sub") cycles to the DoD sequence*

	for i in 1:(length(soc_peak_id)-1)
		push!(DoD_seq, abs( soc[soc_peak_id[i+1]] - soc[soc_peak_id[i]] )) #DoD
		push!(mean_Soc_seq, (soc[soc_peak_id[i+1]] + soc[soc_peak_id[i]]) /2 ) #SoC mean
		push!(delta_t_seq, 0)#soc_peak_id[i+1] - soc_peak_id[i])
	end



	for i in 1:length(DoD_seq)
		Sum_fd += compute_fd(liion.Electro_chimique_params,  DoD_seq[i] , liion.temperature, mean_Soc_seq[i], delta_t_seq[i])
	end

	L = 1 - ( liion.Electro_chimique_params.alpha_sei * exp(-liion.Electro_chimique_params.beta_sei * Sum_fd) ) - ( (1 - liion.Electro_chimique_params.alpha_sei) * exp(-Sum_fd ))

	#SOH = 1- L
	return 1 - L, Sum_fd

end

function S_delta(params::Electro_chimique_params, DoD::Float64)
	return params.k_delta1 * (DoD ^ params.k_delta2)
end

function S_T(params::Electro_chimique_params, T::Float64)
	return exp(params.k_T*(T-params.T_ref) * (params.T_ref/T))
end

function S_sigma(params::Electro_chimique_params, mean_SoC::Float64)
	return exp(params.k_sigma * (mean_SoC  - params.sigma_ref))
end

function S_t(params::Electro_chimique_params, t::Int64)
	return params.k_t * t * 3600 #hours to second
end

function compute_fd(params::Electro_chimique_params, DoD::Float64, T::Float64, mean_SoC::Float64, t::Int64)
	return (0.5 * S_delta(params, DoD) + S_t(params, t)) * S_sigma(params, mean_SoC) * S_T(params, T)
end

 ### Investment dynamic
 function compute_investment_dynamics!(y::Int64, s::Int64, liion::Liion_electro_chimique, decision::Union{Float64, Int64})
     liion.Erated[y+1,s], liion.soc[1,y+1,s], liion.soh[1,y+1,s], liion.voltage[1,y+1,s], liion.Sum_fd[s] = compute_investment_dynamics(liion, (Erated = liion.Erated[y,s], soc = liion.soc[end,y,s], soh = liion.soh[end,y,s], voltage = liion.voltage[end,y,s], Sum_fd = liion.Sum_fd[s]), decision)
 end

 function compute_investment_dynamics(liion::Liion_electro_chimique, state::NamedTuple{(:Erated, :soc, :soh, :voltage, :Sum_fd), Tuple{Float64, Float64, Float64, Float64, Float64}}, decision::Union{Float64, Int64})
     if decision > 1e-2
         Erated_next = decision
         soc_next = liion.soc_ini
         soh_next =  1.
		 Sum_fd_next = 0.
		 #battery replacement
     else
         Erated_next = state.Erated
         soc_next = state.soc
         soh_next = state.soh
		 Sum_fd_next = state.Sum_fd
     end

	 if state.voltage != 0
	 	voltage_next = state.voltage
	 else
		 voltage_next = 3.7
		#error(state.voltage ," is not a value suited for voltage, look like there is a probleme with the computation of this field ")
	 end

     return Erated_next, soc_next, soh_next, voltage_next, Sum_fd_next
 end


#for rule based #TODO change computataion for soc and soh
 function compute_operation_dynamics(liion::Liion_electro_chimique, state::NamedTuple{(:Erated, :soc, :soh), Tuple{Float64, Float64, Float64}}, decision::Float64, Δh::Int64)
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
 	SoH_threshold::Float64, couplage::NamedTuple{(:E,:R), Tuple{Bool,Bool}}, soc_model::String, Erated_ini::Float64, soc_ini::Float64,
 	soh_ini::Float64, update_by_year::Int64, artificial_soc_profil::Array{Float64,2})

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

 	if 8760%update_by_year != 0
 		error("update_by_year value : ", update_by_year, ", do not divide 8760 (365j * 24h). please chose a proper divider")
 		validation = false
 	end

 	return validation
 end
