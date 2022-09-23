abstract type AbstractLiion <: AbstractStorage  end

mutable struct Tremblay_dessaint_params
    E0::Float64
    A::Float64
    B::Float64
    K::Float64
    R::Float64
    Npara::Int64
    Nserie::Int64
	I_max::Float64

    Tremblay_dessaint_params(; E0 = 3.41,#3.366,
                        A = 0.63, #0.26422,
                        B = 4.13e-2,#26.5487,
                        K = 2.44e-4,# 0.0076,
                        R = 3.53e-3,#0.01,
                        Npara = 1,
                        Nserie = 1,
						I_max = 10000.) = new(E0, A, B, K, R, Npara, Nserie, I_max)
end



#TODO generalize every function realtive to the soc

########################
##### SOC models #######
########################
function compute_operation_soc_linear(liion::AbstractLiion, state::NamedTuple{(:Erated, :soc, :soh), Tuple{Float64, Float64, Float64}}, decision::Float64, Δh::Int64)


	η_ini = 0.95   #Fixed (dis)charging efficiency for both BES and EV (0.98)     dans la nomenclature

	if liion.couplage.E
	 Erated = state.Erated * state.soh
	else
	 Erated = state.Erated
	end

	if liion.couplage.R
		η = η_ini - ((1-state.soh)/12)   #(15) simplifié
	else
		η = η_ini
	end

	#power_dch = max(min(decision, liion.α_p_dch * Erated, state.soh * state.Erated / Δh, η * (state.soc * (1. - liion.η_self * Δh) - liion.α_soc_min) * Erated / Δh), 0.)
	#power_ch = min(max(decision, -liion.α_p_ch * Erated, -state.soh * state.Erated / Δh, (state.soc * (1. - liion.η_self * Δh) - liion.α_soc_max) * Erated / Δh / η), 0.)
	#without self discharge
	power_dch, power_ch = get_power_flow(liion, state, decision, Δh)

	#return state.soc * (1. - liion.η_self * Δh) - (power_ch * η + power_dch / η) * Δh / Erated, power_ch + power_dch
	return state.soc - (power_ch * η + power_dch / η) * Δh / Erated, power_ch + power_dch
 end


 function get_power_flow(liion::AbstractLiion, state::NamedTuple{(:Erated, :soc, :soh), Tuple{Float64, Float64, Float64}}, decision::Float64, Δh::Int64)
	 if liion.couplage.E
 	 Erated = state.Erated * state.soh
 	else
 	 Erated = state.Erated
 	end

	η_ini = 0.95   #Fixed (dis)charging efficiency for both BES and EV (0.98)     dans la nomenclature

	if liion.couplage.R
		η = η_ini - ((1-state.soh)/12)   #(15) simplifié
	else
		η = η_ini
	end

	power_dch = max(min(decision, liion.α_p_dch * Erated, state.soh * state.Erated / Δh, η * (state.soc - liion.α_soc_min) * Erated / Δh), 0.)
 	power_ch = min(max(decision, -liion.α_p_ch * Erated, -state.soh * state.Erated / Δh, (state.soc - liion.α_soc_max) * Erated / Δh / η), 0.)

	return power_dch, power_ch
end


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

	 Qnom = state.Erated / (3.7 * liion.tremblay_dessaint_params.Npara * liion.tremblay_dessaint_params.Nserie) # TODO 3.7 en V pour convertir Erated de kWh en Ah
	 if liion.couplage.E
     	Qnom *= state.soh #TODO faire verifier
	 end


     I = decision / (V * liion.tremblay_dessaint_params.Npara * liion.tremblay_dessaint_params.Nserie)
     it = (1-state.soc)*Qnom

     # calcul du soc suivant et volt suivant

     it_next = it + I * Δh
     soc_next = 1 - it_next/Qnom



	 soc_next = min(max(soc_next,soc_min),soc_max)
	 it_next = (1-soc_next)*Qnom
	 I = (it_next - it)/Δh


	 decision = I * V
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




#Optimal Sizing and Control of a PV-EV-BES Charging System Including Primary Frequency Control and Component Degradation
#Wiljan Vermeer et al.
#With this soc the efficiency is based on battery state of health
 function compute_operation_soc_Vermeer(liion::AbstractLiion, state::NamedTuple{(:Erated, :soc, :soh), Tuple{Float64, Float64, Float64}}, decision, Δh::Int64)

	η_ini = 0.98   #Fixed (dis)charging efficiency for both BES and EV (0.98)     dans la nomenclature
	η = η_ini - ((1-state.soh)/12)   #(15) simplifié

	power_dch = max(min(decision, liion.α_p_dch * state.Erated, state.soh *  state.Erated / Δh, η * (liion.α_soc_max - state.soc) * state.Erated * state.soh / Δh), 0.)
 	power_ch = min(max(decision, -liion.α_p_ch * state.Erated, -state.soh *  state.Erated / Δh, (liion.α_soc_min - state.soc) * state.Erated * state.soh / Δh / η), 0.)

	P = ( power_dch / η) + (power_ch * η) #(13)
	E_lim = state.Erated * state.soh   #Definition : Maximum battery capacity at time t, based on degradation
	E = E_lim * state.soc # based on (31)
	new_E = E + P * Δh #(34)

	return new_E/E_lim , P# Get the SoC to keep coherence with the entire code.
  end

function compute_operation_soc_artificial(liion::AbstractLiion, profil::Array{Float64,2},  y::Int64, h::Int64)
	return profil[y,h]
end



function get_soc_peaks(soc::Vector{Float64})
	soc_peak = Float64[] #soc_peak is the sequence of values for the state of charges peaks
	soc_peak_id = Int64[] #soc_peak is the sequence of values for the state of charges peaks

	#add first
	push!(soc_peak, soc[1])
	push!(soc_peak_id, 1)

	# Extract all peak value of soc
	for i in 2:(length(soc)-1)
		sign1 = soc[i] - soc[i-1] > 0 #true = positiv, false = negativ
		sign2 = soc[i+1] - soc[i] > 0

		if sign1 != sign2 # different sign mean change of trend (increasing of decreasing) so it's a peak
			push!(soc_peak, soc[i])
			push!(soc_peak_id, i)
		end
	end



	#add last
	push!(soc_peak, soc[length(soc)])
	push!(soc_peak_id, length(soc))

	return soc_peak, soc_peak_id
end
