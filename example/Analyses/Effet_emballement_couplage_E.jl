
#Effet d'emballement du soh avec le couplage E
#Comment ça se dégrade différement et pourquoi.

#On observe plusieurs phénomènes intéréssants.

#Déjà avec Rainflow lorsqu'on augmente la taille de la batterie on obtient une dégradation plus rapide de celle ci (ce qui est contre intuitif) mais cela est du à une utilisaation plus importante de celle ci
#Cela est vrai jusqu'a un certain seuil où la batterie sera solicitée pour les mêmes raisons mais dès lors les DoD seront de plus en plus petit et alors la batterie se dégradera de moins en moins vite avec l'augmentation de sa taille

#On observe aussi l'emballement  du vieillissement pour des batteries de grandes tailles qui ont originnellement une capacité suffisemment importante pour les cycles ne soient généralement pas de borne à borne.
#Le vieillissement augmente la profondeur des cycles lorsqu'il y a de la marge jusqu'a un certain points où les cycles sont maximum et la batterie est de moins en moins utilisé, la plongé arrête alors et on retrouve une dégradation quasi linéaire.

#hypothèse : lorsque la profondeur des cycles croit elle le fait selon une loi bruité en x^2 du temps

#Les différences de modèle énergétique sont notables au niveau du rendement, tremblay_dessaint renvoi un rendement plus faibnle ce qui a tendance à préserver la batterie du vieillissement avec des cycles sensiblement aussi nombreux mais moins profond.

#Question supplémentaire : qu'en est-il de cet impact sur l'autonomie

# Load packages
include("..\\..\\src\\Genesys.jl")

using Main.Genesys
using JLD, Dates, Seaborn
using CSV, DataFrames
using Statistics

pygui(true)

# Parameters of the simulation
const nh, ny, ns = 8760, 8, 1

# Load input data
data = load(joinpath("example","data","ausgrid_5_twostage.jld"))

#https://www.researchgate.net/figure/Battery-cycle-to-failure-versus-the-DOD-gap-34_fig5_276136726  for data1
# Modeling of Lithium-Ion Battery Degradation for Cell Life Assessment, Bolun Xu et al. for data 2
fatigue_data = DataFrame(CSV.File("example\\data\\fatigue_data_NMC.csv", delim = ",", header = [Symbol("DoD"),Symbol("cycle")], types=Dict(:DoD=>Float64, :cycle=>Float64)))



mg = Microgrid(parameters = GlobalParameters(nh, ny, ns, renewable_share = 1.))

#Create microgrid
add!(mg, Demand(carrier = Electricity()),
                Solar(),
                Liion_electro_chimique(soc_model = "linear"),
                Grid(carrier = Electricity()))

ω_d = Scenarios(mg, data["ω_optim"]; same_year=true, seed = 25:25)


using Sobol

PV_range = [5,50]
BAT_range = [10,80]

P_sous = 20

Nb_sizing = 512

lb = [PV_range[1], BAT_range[1], P_sous]
ub = [PV_range[2], BAT_range[2], P_sous]
seq = SobolSeq(lb, ub)
sizings = transpose(Base.reduce(hcat, next!(seq) for i = 1:Nb_sizing))
sizings_tuple = [(pv = sizings[i,1], bat = sizings[i,2], P_sub = sizings[i,3]) for i in 1:Nb_sizing]


Seaborn.scatter(sizings[:,1], sizings[:,2])

function get_DoD_seq(soc_profil::Vector{Float64})

    soc_peak = Float64[] #soc_peak is the sequence of values for the state of charges peaks

    #add first
    push!(soc_peak, soc_profil[1])
    # Extract all peak value of soc
    for i in 2:(length(soc_profil)-1)
        sign1 = soc_profil[i] - soc_profil[i-1] > 0 #true = positiv, false = negativ
        sign2 = soc_profil[i+1] - soc_profil[i] > 0

        if sign1 != sign2 # different sign mean change of trend (increasing of decreasing) so it's a peak
            push!(soc_peak, soc_profil[i])
        end
    end

    #add last
    push!(soc_peak, soc_profil[length(soc_profil)])
    peak_copy = copy(soc_peak)
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
            deleteat!(soc_peak, i+1)
            deleteat!(soc_peak, i+2)
        else #else use the following point sequence
            i = i+1
        end
    end

    #Then add the englobing (those who make the other cycles "sub") cycles to the DoD sequence*
    for i in 1:(length(soc_peak)-1)
        push!(DoD_seq, abs(soc_peak[i+1]-soc_peak[i]))
    end


	#currently neglect cycle under 1%
	deleteat!(DoD_seq, findall(<(1e-2), DoD_seq))


    return DoD_seq
end

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


using Pandas
using PyPlot
using PyCall

@pyimport numpy

function Φ(DoD::Float64, fatigue_data)
	 index = findfirst(>=(DoD), fatigue_data.DoD)

	 return fatigue_data.cycle[index]
end

#Initiate dataframe
df = Pandas.DataFrame(Dict(:SoH=>[], :sizing_num=>[], :DoDs=>[],
 :SoC_moyen=>[], :Cycles=>[], :Fatigue=>[], :SoC_mean=>[],
  :energie=>[], :RES=>[], :couplage=>[], :energie_bought=>[]))


for couplage in [(E=false,R=false), (E=true,R=false), (E=false,R=true), (E=true,R=true)]


	fatigues = zeros(Nb_sizing, ny-1)
	DoD_mean = zeros(Nb_sizing, ny-1)
	SoC_mean = zeros(Nb_sizing, ny-1)
	N_cycle = zeros(Nb_sizing, ny-1)
	energie = zeros(Nb_sizing, ny-1)
	RES = zeros(Nb_sizing, ny-1)
	energie_bought = zeros(Nb_sizing, ny-1)

	for (k,sizes) in enumerate(sizings_tuple)


		param = [sizes.pv, sizes.bat, sizes.P_sub]
		controller = RBC(options = RBCOptions(policy_selection = 2 ))
		designer = Manual(generations = [param[1]], storages = [param[2]], subscribed_power = [param[3]])

		battery = Liion_vermeer(update_by_year = 12, soc_model = "linear", couplage = couplage, SoH_threshold = 0.3)

		mg = Microgrid(parameters = GlobalParameters(nh, ny, ns, renewable_share = 1.))

		#Create microgrid
		add!(mg, Demand(carrier = Electricity()),
		                Solar(),
		                battery,
		                Grid(carrier = Electricity()))



		#generate scenarios  #Too many for the moment by the way, the need is different for this study

		designer = initialize_designer!(mg, designer, ω_d)
		controller = initialize_controller!(mg, controller, ω_d)

		simulate!(mg, controller, designer, ω_d, options = Genesys.Options(mode = "serial", firstyear = false))
		metrics = Metrics(mg, designer)
		#plot_operation(mg, y=2:ny, xdisplay = "years", ci = true)





		# years = []
		# scenarios = []
		# DoD_mean_by_years = []
		# SoC_mean_by_years = []
		# N_cycle_by_years = []
		# fatigue_by_year = []
		# energie_by_year = []
		# RES_by_year = []
		# couplages = []
		# energie_bought_by_year = []
		#
		# DoD_for_rempl = zeros(ny-1)
		# SoC_for_rempl = zeros(ny-1)
		# N_cycle_for_rempl = zeros(ny-1)
		# fatigue_dor_rempl = zeros(ny-1)
		# energie_for_rempl = zeros(ny-1)
		# RES_for_rempl = zeros(ny-1)
		# energie_bought_for_rempl = zeros(ny-1)

		for y in 2:ny


			DoD_seq = get_DoD_seq(vec(mg.storages[1].soc[:,y,1]))

	 		fatigue = 0
 		   	for i in 1:length(DoD_seq)
 			   	fatigue += 1/(2*Φ(DoD_seq[i], fatigue_data)) #Compute fatigue with phy function applied to all the half cycles DoD factor 2 refer to half cycles
 		   	end

	 		fatigues[k,y-1] = fatigue
	 		DoD_mean[k,y-1] = mean(DoD_seq)
	 		N_cycle[k,y-1] = length(DoD_seq)


			energie[k,y-1] = 0


			soh_coef = ones(nh+1)
			if couplage.E
				soh_coef .= mg.storages[1].soh[:,y,1]
			end

			for t in 1:(length(vec(mg.storages[1].soc[:,y,1]))-1)
				energie[k,y-1] += abs(mg.storages[1].soc[t,y,1] - mg.storages[1].soc[t+1,y,1]) * sizes[2] * soh_coef[t]
			end
			RES[k,y-1] = metrics.renewable_share[y,1] * 100

			SoC_mean[k,y-1] = mean(mg.storages[1].soc[:,y,1])

			energie_bought[k,y-1] = sum(max.(0.,  mg.grids[1].carrier.power[:,y,1]))


		end

		# years = vcat(years,repeat([y],ns))
		#
		# 	couplages = vcat(couplages,repeat([couplage],ns))
		#
		# 	scenarios = vcat(scenarios, [a for a in 1:ns])
		#
		#  	DoD_mean_by_years = vcat(DoD_mean_by_years, DoD_mean)
		# 	SoC_mean_by_years = vcat(SoC_mean_by_years, SoC_mean)
		# 	N_cycle_by_years = vcat(N_cycle_by_years, N_cycle)
		# 	fatigue_by_year = vcat(fatigue_by_year, fatigues)
		# 	energie_by_year = vcat(energie_by_year, energie)
		# 	RES_by_year = vcat(RES_by_year, RES)
		# 	energie_bought_by_year = vcat(energie_bought_by_year, energie_bought)
		#
		# 	DoD_for_rempl[y-1] = mean(DoD_mean)
		# 	SoC_for_rempl[y-1] = mean(SoC_mean)
		# 	N_cycle_for_rempl[y-1] = mean(N_cycle)
		# 	fatigue_dor_rempl[y-1] = mean(fatigues)
		# 	energie_for_rempl[y-1] = mean(energie)
		# 	RES_for_rempl[y-1] = mean(RES)
		# 	energie_bought_for_rempl[y-1] = mean(energie_bought)





	end


	SoH = repeat(mg.storages[1].soh[1,2:ny,1], inner=Nb_sizing)
	sizing_num = repeat(1:Nb_sizing,ny-1)
	couplage = repeat([string(couplage)], Nb_sizing*(ny-1) )


	DoD_mean_percent = 100 .* (DoD_mean .- DoD_mean[:,1]) ./ DoD_mean[:,1]
	SoC_mean_percent = 100 .* (SoC_mean .- SoC_mean[:,1]) ./ SoC_mean[:,1]
	energie_percent =  100 .* (energie .- energie[:,1]) ./ energie[:,1]
	energie_bought_percent = 100 .* (energie_bought .- energie_bought[:,1]) ./ energie_bought[:,1]

	df_tmp = Pandas.DataFrame(Dict(:SoH=>SoH, :sizing_num=>sizing_num, :DoDs=>vec(DoD_mean_percent),
	 :SoC_moyen=>vec(SoC_mean_percent), :Cycles=>vec(N_cycle), :Fatigue=>vec(fatigues),
	  :energie=>vec(energie_percent), :SoC_mean=>vec(SoC_mean), :RES=>vec(RES), :couplage=>couplage, :energie_bought=>vec(energie_bought_percent)))


	if isempty(df)
     df = df.pyo[:append](df_tmp, ignore_index = true )
  	else
     df = df[:append](df_tmp, ignore_index = true )
  	end



end



#
# remplacement = (mg.storages[1].soh[end,1:ny,1:ns] .< mg.storages[1].SoH_threshold)
# remplacement  = [sum(remplacement[a,:]) * 200/ns for a in 2:(ny)]
# df_rempl = Pandas.DataFrame(Dict(:Années=>1:(ny-1), :DoDs=>DoD_for_rempl, :SoC_moyen=>SoC_for_rempl, :fatigue_by_year => fatigue_dor_rempl, :N_cycle_by_years => N_cycle_for_rempl,:energie=>energie_for_rempl, :remplacement=>remplacement, :RES=>RES_for_rempl, :energie_bought=>energie_bought_for_rempl))
# df2_rempl = query(df_rempl, :(remplacement>0))# or query(df3, "income>85")
# df2_rempl[:Années] = df2_rempl[:Années].+1
#
# x_ticks_lab = [string(n) for n in 1:(ny-1)]
# x_ticks = [n for n in 2:(ny)]


axeslabelsize = 31
legendfontsize = 31
legendtitle_fontsize = 31
yticklabelsize = 34
xticklabelsize = 34

figure(string("Mean DoD evolution over years"))
Seaborn.set_context("paper", rc=Dict("axes.labelsize"=>axeslabelsize, "legend.fontsize"=> legendfontsize,  "legend.title_fontsize"=> legendtitle_fontsize,
"ytick.labelsize" => yticklabelsize, "xtick.labelsize" => xticklabelsize))
plotDoD = Seaborn.lineplot(data=df, x="SoH", y="DoDs", hue="couplage", linewidth = 2, markers=true, dashes=false, ci=100)
xlabel("SoH")
ylabel("Annual average DoD of cycles (%)")
legend(["None","E", "R", "ER"], title="Coupling")

plotDoD.invert_xaxis()
#
# Seaborn.scatter(data=df2_rempl, x="Années", y="DoDs", marker="^", color = "red", s=:remplacement)
# #legend([string("Erated = ", 10 , " kWh"), string("Erated = ", 50 , " kWh"), string("Erated = ", 100 , " kWh"),string("Erated = ", 150 , " kWh"), "CI 95/5" , string("Rempl. Batterie")], loc = 3)
# legend([string((E=false,R=false)), string((E=true,R=false)), string((E=false,R=true)), string((E=true,R=true)), "CI 95/5" , string("Rempl. Batterie")], loc = 3)
# plotDoD.set_xticks(x_ticks)
# plotDoD.set_xticklabels(x_ticks_lab)

figure(string("Evolution of the average SoC of cycles (%)"))
Seaborn.set_context("paper", rc=Dict("axes.labelsize"=>axeslabelsize, "legend.fontsize"=> legendfontsize,  "legend.title_fontsize"=> legendtitle_fontsize,
"ytick.labelsize" => yticklabelsize, "xtick.labelsize" => xticklabelsize))
plotSoC = Seaborn.lineplot(data=df, x="SoH", y="SoC_moyen", hue="couplage",
 linewidth = 2, markers=true, dashes=false, ci=100)
xlabel("SoH")
ylabel("Annual average SoC of cycles (%)")
legend(["None","E", "R", "ER"], title="Coupling", loc="upper left")



plotSoC.invert_xaxis()
#
# Seaborn.scatter(data=df2_rempl, x="Années", y="SoC_moyen", marker="^", color = "red", s=:remplacement)
# #legend([string("Erated = ", 10 , " kWh"), string("Erated = ", 50 , " kWh"), string("Erated = ", 100 , " kWh"),string("Erated = ", 150 , " kWh"), "CI 95/5" , string("Rempl. Batterie")], loc = 3)
# legend([string((E=false,R=false)), string((E=true,R=false)), string((E=false,R=true)), string((E=true,R=true)), "CI 95/5" , string("Rempl. Batterie")], loc = 3)
# plotDoD.set_xticks(x_ticks)
# plotDoD.set_xticklabels(x_ticks_lab)

figure(string("Energy exchanged over years"))
Seaborn.set_context("paper", rc=Dict("axes.labelsize"=>axeslabelsize, "legend.fontsize"=> legendfontsize,  "legend.title_fontsize"=> legendtitle_fontsize,
"ytick.labelsize" => yticklabelsize, "xtick.labelsize" => xticklabelsize))
plotEE = Seaborn.lineplot(data=df, x="SoH", y="energie", hue="couplage", linewidth = 2, ci = 100)
xlabel("SoH")
ylabel("Annual exchanged energy (%)")
legend(["None","E", "R", "ER"], title="Coupling")

plotEE.invert_xaxis()


# Seaborn.scatter(data=df2_rempl, x="Années", y="energie", marker="^", color = "red", s=:remplacement)
# legend([string((E=false,R=false)), string((E=true,R=false)), string((E=false,R=true)), string((E=true,R=true)), "CI 95/5" , string("Rempl. Batterie")], loc = 3)
# plotEE.set_xticks(x_ticks)
# plotEE.set_xticklabels(x_ticks_lab)



figure(string("Energy bought over years"))
Seaborn.set_context("paper", rc=Dict("axes.labelsize"=>axeslabelsize, "legend.fontsize"=> legendfontsize,  "legend.title_fontsize"=> legendtitle_fontsize,
"ytick.labelsize" => yticklabelsize, "xtick.labelsize" => xticklabelsize))
plotEB = Seaborn.lineplot(data=df, x="SoH", y="energie_bought", hue="couplage", linewidth = 2, ci = 100)
xlabel("SoH")
ylabel("Annual purchased energy (%)")
legend(["None","E", "R", "ER"], title="Coupling")


plotEB.invert_xaxis()





figure(string("RES"))
Seaborn.set_context("paper", rc=Dict("axes.labelsize"=>30, "legend.fontsize"=> 30,  "legend.title_fontsize"=> 30,
"ytick.labelsize" => 25, "xtick.labelsize" => 25))
plotEB = Seaborn.lineplot(data=df, x="SoH", y="RES", hue="couplage", linewidth = 2, ci = 95)
xlabel("SoH")
ylabel("RES")
legend(title="Couplage")

plotEB.invert_xaxis()


# Seaborn.scatter(data=df2_rempl, x="Années", y="energie_bought", marker="^", color = "red", s=:remplacement)
# legend([string((E=false,R=false)), string((E=true,R=false)), string((E=false,R=true)), string((E=true,R=true)), "CI 95/5" , string("Rempl. Batterie")], loc = 3)
# plotEE.set_xticks(x_ticks)
# plotEE.set_xticklabels(x_ticks_lab)



figure("Number of cycle evolution over year")
plotCycle = Seaborn.lineplot(data=df, x="Years", y="Cycles", linewidth = 2, ci = 95)
xlabel("Années")
ylabel("Cycles/an")
Seaborn.scatter(data=df2_rempl, x="Années", y="N_cycle_by_years", marker="^", color = "red", s=:remplacement)
legend([string("Erated = ", 10 , " kWh"), string("Erated = ", 50 , " kWh"), string("Erated = ", 100 , " kWh"),string("Erated = ", 150 , " kWh"), "CI 95/5" , string("Rempl. Batterie")], loc = 3)
plotCycle.set_xticks(x_ticks)
plotCycle.set_xticklabels(x_ticks_lab)

figure("Fatigue over years")
plotCycle = Seaborn.lineplot(data=df, x="Years", y="Fatigue", linewidth = 2, ci = 95)
xlabel("Années")
ylabel("Fatigue (perte de SoH)")
Seaborn.scatter(data=df2_rempl, x="Années", y="fatigue_by_year", marker="^", color = "red", s=:remplacement)
legend([string("Erated = ", 10 , " kWh"), string("Erated = ", 50 , " kWh"), string("Erated = ", 100 , " kWh"),string("Erated = ", 150 , " kWh"), "CI 95/5" , string("Rempl. Batterie")], loc = 3)
plotCycle.set_xticks(x_ticks)
plotCycle.set_xticklabels(x_ticks_lab)

figure("Energy exchanged over years")
plotEE = Seaborn.lineplot(data=df, x="Years", y="energie", linewidth = 2, ci = 95)
xlabel("Années")
ylabel("Energie echangée (kWh)")
Seaborn.scatter(data=df2_rempl, x="Années", y="energie", marker="^", color = "red", s=:remplacement)
legend([string("Erated = ", 10 , " kWh"), string("Erated = ", 50 , " kWh"), string("Erated = ", 100 , " kWh"),string("Erated = ", 150 , " kWh"), "CI 95/5" , string("Rempl. Batterie")], loc = 3)
plotEE.set_xticks(x_ticks)
plotEE.set_xticklabels(x_ticks_lab)


figure("Autonomie over years")
plotRES = Seaborn.lineplot(data=df, x="Years", y="RES", linewidth = 2, ci = 95)
xlabel("Années")
ylabel("Autonomie energétique (%)")
Seaborn.scatter(data=df2_rempl, x="Années", y="RES", marker="^", color = "red", s=:remplacement)
legend([string("Erated = ", 10 , " kWh"), string("Erated = ", 50 , " kWh"), string("Erated = ", 100 , " kWh"),string("Erated = ", 150 , " kWh"), "CI 95/5" , string("Rempl. Batterie")], loc = 3)
plotRES.set_xticks(x_ticks)
plotRES.set_xticklabels(x_ticks_lab)


# PyPlot.plot(N_cycle_by_years)
# xlabel("Years")
# ylabel("Number of cycle")
# legend([string("bat = ", 10 , " kWh"),string("bat = ", 50 , " kWh"),string("bat = ", 150 , " kWh"),string("bat = ", 300 , " kWh"), string("bat = ", 300 , " kWh, threshhold 0.8")], loc = 3)
# figure("DoD over years")
# PyPlot.plot(DoD_mean_by_years)
# legend([string("bat = ", 10 , " kWh"),string("bat = ", 50 , " kWh"),string("bat = ", 150 , " kWh"),string("bat = ", 300 , " kWh"), string("bat = ", 300 , " kWh, threshhold 0.8")], loc = 3)
# figure("fatigue over years")
# PyPlot.plot(fatigue_by_year)
# legend([string("bat = ", 10 , " kWh"),string("bat = ", 50 , " kWh"),string("bat = ", 150 , " kWh"),string("bat = ", 300 , " kWh"), string("bat = ", 300 , " kWh, threshhold 0.8")], loc = 3)
#



print( get_battery_lifetime(mg.storages[1]))


#df2 = Pandas.DataFrame(Dict(:years=>(2:ny), :Fatigue_mean=>fatigue_by_year))


#lmplot(x="years", y="Fatigue_mean", data=df2[1:(get_battery_lifetime(mg.storages[1])-2)], order =2);
#lmplot(x="years", y="DoD_med", data=df, order =2);


Seaborn.plot(900:1000, 200 ./(1000:-2:800))
