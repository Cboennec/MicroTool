#=
    Li-ion battery modelling
 =#


 soc_model_names = ["tremblay_dessaint", "linear", "vermeer", "artificial"]

#rainflow ref : Optimal Battery Control Under Cycle Aging Mechanisms in Pay for Performance Settings
# Yuanyuan Shi, Bolun Xu, Yushi Tan, Daniel Kirschen, Baosen Zhang
mutable struct Liion_rainflow <: AbstractLiion
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

	#Model dynamics
	soc_model::String #model name
	calendar::Bool

	# Initial conditions
	Erated_ini::Float64  # capacité de la batterie en Wh
	soc_ini::Float64 # first state of charge for the begining of simulation
	soh_ini::Float64 # first state of health for the begining of simulation

	#input specific params
	update_by_year::Int64 # Rainflow Soh computation by year

	fatigue_data::DataFrame #2 column data frame (DoD, ncycle) giving the number of cycle that the battery can do for a given DoD

	artificial_soc_profil::Array{Float64,2} # soc profil for study purpose

	# Variables
	Erated::AbstractArray{Float64,2} # Battery capacity
	carrier::Electricity #Type of energy
	soc::AbstractArray{Float64,3} #3 dim matrix (h,y,s) containing the state of charge [0-1]
	soh::AbstractArray{Float64,3} #3 dim matrix (h,y,s) containing the state of health [0-1]
	#tremblay_dessaint
	voltage ::Array{Float64,3} # evolution de la tension
	current ::Array{Float64,3} # evolution du courant ?
	tremblay_dessaint_params::Tremblay_dessaint_params # A list of params only use for tremblay dessaint


	#inner specific params

	# Eco
	cost::AbstractArray{Float64,2}

	# Inner constructor
	Liion_rainflow(; α_p_ch = 1.5,
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
		soc_model = "linear",
		calendar = false,
		Erated_ini = 1e-6,
		soc_ini = 0.5,
		soh_ini = 1.,
		update_by_year = 12,
		fatigue_data = DataFrames.DataFrame(CSV.File("example\\data\\fatigue_data2.csv.csv", delim = ";", header = [Symbol("DoD"),Symbol("cycle")], types=Dict(:DoD=>Float64, :cycle=>Float64))),
		artificial_soc_profil = zeros(8760,1)
		) =  verification_liion_params(α_p_ch, α_p_dch, η_ch, η_dch, η_self, α_soc_min, α_soc_max, lifetime, nCycle, bounds,
			SoH_threshold, couplage, soc_model, calendar, Erated_ini, soc_ini, soh_ini, update_by_year, artificial_soc_profil) ?
			new(α_p_ch, α_p_dch, η_ch, η_dch, η_self, α_soc_min, α_soc_max, lifetime, nCycle, bounds,
			SoH_threshold, couplage, soc_model, calendar, Erated_ini, soc_ini, soh_ini, update_by_year, fatigue_data, artificial_soc_profil) : nothing
end

### Preallocation
function preallocate!(liion::Liion_rainflow, nh::Int64, ny::Int64, ns::Int64)
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
	#liion.fatigue_data[:,"cycle"] *= 1
	liion.tremblay_dessaint_params = Tremblay_dessaint_params()

	return liion
end

### Operation dynamic
function compute_operation_dynamics!(h::Int64, y::Int64, s::Int64, liion::Liion_rainflow, decision::Float64, Δh::Int64)

	if liion.soc_model == "tremblay_dessaint"
		liion.soc[h+1,y,s], liion.voltage[h+1,y,s], liion.carrier.power[h,y,s], liion.current[h,y,s] = compute_operation_soc_tremblay_dessaint(liion, (Erated = liion.Erated[y,s], soc = liion.soc[h,y,s], soh = liion.soh[h,y,s]),  liion.voltage[h,y,s], decision, Δh)
	elseif  liion.soc_model == "linear"
		liion.soc[h+1,y,s], liion.carrier.power[h,y,s] = compute_soc_dynamics(liion, (Erated = liion.Erated[y,s], soc = liion.soc[h,y,s], soh = liion.soh[h,y,s]), decision, Δh)
	elseif liion.soc_model == "artificial"

	end


#	if liion.soc[h+1,y,s] > 0.805 || liion.soc[h+1,y,s] < 0.195
#		println("soc bound break : ", liion.soc[h+1,y,s] )
#	end

	h_between_update = convert(Int64,floor(8760/liion.update_by_year))

	#SoH computation
    if (h%h_between_update) != 0
		liion.soh[h+1,y,s] = liion.soh[h,y,s]
    else #rainflow computaion
		interval = (h-h_between_update+1):h

		if isnan(liion.soc[interval[2],y,s])
			println("y = : ",y, ";     soc = ", soc)
			println("interval = : ", interval)
		end

		liion.soh[h+1,y,s] = compute_operation_soh_rainflow(liion, (Erated = liion.Erated[y,s], soc = liion.soc[h,y,s], soh = liion.soh[h,y,s]), decision, Δh,  liion.soc[interval,y,s])

		#Calendar part
		if liion.calendar == true
            liion.soh[h+1,y,s] = liion.soh[h+1,y,s] - (1 - exp(- 4.14e-10 * 3600 * h_between_update))
        end
    end

end

function compute_soc_dynamics(liion::Liion_rainflow, state::NamedTuple{(:Erated, :soc, :soh), Tuple{Float64, Float64, Float64}}, decision::Float64, Δh::Int64)

    # SoC dynamic
    return compute_operation_soc_linear(liion, state, decision, Δh)

end






########################
##### SOH models #######
########################

function Φ(DoD::Float64, fatigue_data)
	 index = findfirst(>=(DoD), fatigue_data.DoD)

     return fatigue_data.cycle[index]
end





function compute_operation_soh_rainflow(liion::Liion_rainflow, state::NamedTuple{(:Erated, :soc, :soh), Tuple{Float64, Float64, Float64}}, decision::Float64, Δh::Int64, soc::Vector{Float64})

	#Gather peaks from the soc profil
	soc_peak, _ = get_soc_peaks(soc)


	#Then compute the DoD sequence by extracting the subcycles DoD
	DoD_seq = Float64[] #Sequence of all the charging and decharging half cycles DoDs

	i = 1

	while i+3 <= length(soc_peak)
		#Define your 3 deltas with 4 consecutives points
		delta1 = abs( soc_peak[i+1] - soc_peak[i] )
		delta2 = abs( soc_peak[i+2] - soc_peak[i+1] )
		delta3 = abs( soc_peak[i+3] - soc_peak[i+2] )

		#rainflow sub-cycle criterion
		if delta2 <= delta1 && delta2 <= delta3
			push!(DoD_seq, delta2) #1 half cycle of DoD delta2 +
			push!(DoD_seq, delta2) #1 half cycle
			deleteat!(soc_peak, i+2) #start with the second or you will delete i+1 and i+3
			deleteat!(soc_peak, i+1)
		else #else use the following point sequence
			i = i+1
		end
	end

	#Then add the englobing (those who make the other cycles "sub") cycles to the DoD sequence*
	for i in 1:(length(soc_peak)-1)
		push!(DoD_seq, abs(soc_peak[i+1]-soc_peak[i]))
	end

	#currently neglect cycle under 1%
	#deleteat!(DoD_seq, findall(<(1e-2), DoD_seq))

	if length(DoD_seq) > 0 && isnan(DoD_seq[1])
		println("soc : ", soc)
	end
	fatigue = 0

	for i in 1:length(DoD_seq)
		if DoD_seq[i] > 0.65 && liion.soc_model != "artificial"
			#println(peak_copy)
		end
		fatigue += 1/(2*Φ(DoD_seq[i], liion.fatigue_data) ) #Compute fatigue with phy function applied to all the half cycles DoD factor 2 refer to half cycles
	end

	return state.soh - (fatigue / 5) #/5 car la courbe cycle to failure donne le nombre de cycle jusqu'à 80% SOH
end


 ### Investment dynamic
 function compute_investment_dynamics!(y::Int64, s::Int64, liion::Liion_rainflow, decision::Union{Float64, Int64})
     liion.Erated[y+1,s], liion.soc[1,y+1,s], liion.soh[1,y+1,s], liion.voltage[1,y+1,s] = compute_investment_dynamics(liion, (Erated = liion.Erated[y,s], soc = liion.soc[end,y,s], soh = liion.soh[end,y,s], voltage = liion.voltage[end,y,s]), decision)
 end

 function compute_investment_dynamics(liion::Liion_rainflow, state::NamedTuple{(:Erated, :soc, :soh, :voltage), Tuple{Float64, Float64, Float64, Float64}}, decision::Union{Float64, Int64})
     if decision > 1e-2
         Erated_next = decision
         soc_next = liion.soc_ini
         soh_next =  1.
     else
         Erated_next = state.Erated
         soc_next = state.soc
         soh_next = state.soh
     end

	 #TODO faire verifier
	 if state.voltage != 0
		 voltage_next = state.voltage
	 else
		voltage_next = liion.soc_ini*0.7 + 3.3
	 end

     return Erated_next, soc_next, soh_next, voltage_next
 end


#for rule based policy 1 #TODO change computataion for soc and soh
function compute_operation_dynamics(liion::Liion_rainflow, state::NamedTuple{(:Erated, :soc, :soh), Tuple{Float64, Float64, Float64}}, decision::Float64, Δh::Int64)
	# Control power constraint and correction
	if liion.couplage.E
		Erated = state.Erated * state.soh
	else
		Erated = state.Erated
	end

	power_dch = max(min(decision, liion.α_p_dch * Erated, state.soh * Erated / Δh, liion.η_dch * (state.soc * (1. - liion.η_self * Δh) - liion.α_soc_min) * Erated / Δh), 0.)
	power_ch = min(max(decision, -liion.α_p_ch * state.Erated, -state.soh * Erated / Δh, (state.soc * (1. - liion.η_self * Δh) - liion.α_soc_max) * Erated / Δh / liion.η_ch), 0.)

	#for rule base, doesnt matter
	soc_next = state.soc * (1. - liion.η_self * Δh) - (power_ch * liion.η_ch + power_dch / liion.η_dch) * Δh / state.Erated
	  # SoH dynamic
	soh_next = state.soh - (power_dch - power_ch) * Δh / (2. * liion.nCycle * (liion.α_soc_max - liion.α_soc_min) * state.Erated)
	return soc_next, soh_next, power_dch + power_ch
end




function verification_liion_params(α_p_ch::Float64, α_p_dch::Float64, η_ch::Float64, η_dch::Float64, η_self::Float64,
	α_soc_min::Float64, α_soc_max::Float64, lifetime::Int64, nCycle::Float64, bounds::NamedTuple{(:lb, :ub), Tuple{Float64, Float64}},
	SoH_threshold::Float64, couplage::NamedTuple{(:E,:R), Tuple{Bool,Bool}}, soc_model::String, calendar::Bool, Erated_ini::Float64, soc_ini::Float64,
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
