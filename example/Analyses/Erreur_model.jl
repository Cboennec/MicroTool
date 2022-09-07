# Cette étude cherche à mesurer l'erreur dans les métriques (RES,NPV) vis à vis d'une référence
# La référence utilisée est notre modèle considéré le plus précis. Ce dernier étant le modèle Electrochimique - tremblay_dessaint
# Cette erreur est mesurer en simulant ns scénario sur n_sizing dimensionnement générer par séquence de sobol dans des bornes fixés.
# La puissance souscrite est fixée.
# l'erreur est mesuré avec (value-ref/ref)





include("..\\..\\src\\Genesys.jl")


using Main.Genesys
using JLD, Dates, Seaborn
using CSV
using Pandas, DataFrames
using Profile
using Sobol

pygui(true)

const nh, ny, ns = 8760, 21, 20

data = load(joinpath("example","data","ausgrid_5_twostage.jld"))

PV_range = [5,50]
BAT_range = [10,80]

P_sous = 20

Nb_sizing = 50
couplages_names = [" ", "E", "R", "ER"]
couplages = [(E=false,R=false), (E=true,R=false), (E=false,R=true), (E=true,R=true)]
SoC_model_display = ["Linéaire", "Tremblay-Dessaint"]
SoH_model_display = ["Énergie échangée", "Rainflow", "Éléctro-chimique", "Fixed"]
SoC_model_display_EN = ["Linear", "Tremblay-Dessaint"]
SoH_model_display_EN = ["Energy Exchanged", "Rainflow", "Electro-Chemical", "Fixed"]
SoC_model = ["linear", "tremblay_dessaint"]
SoH_model = ["energie echangee", "rainflow", "electro-chimique", "fixed"]

metriques_names_list = ["Cost", "RES", "NPV*", "Discounted opex", "Discounted capex", "Discounted salvage", "NPV"]

   value = []
   pv = []
   bat = []
   Couplage = []
   SOH = []
   SOC = []
   metriques_names = []

   lb = [PV_range[1], BAT_range[1]]
   ub = [PV_range[2], BAT_range[2]]
   seq = SobolSeq(lb, ub)
   sizings = transpose(Base.reduce(hcat, next!(seq) for i = 1:Nb_sizing))

   #To vizualise param distribution

   #subplot(111, aspect="equal")
   #plot(sizings[:,1], sizings[:,2], "r.")

   fatigue_data = DataFrames.DataFrame(CSV.File("example\\data\\fatigue_data_NMC.csv", delim = ",", header = [Symbol("DoD"),Symbol("cycle")], types=Dict(:DoD=>Float64, :cycle=>Float64)))

   mg = Microgrid(parameters = GlobalParameters(nh, ny, ns, renewable_share = 1.))
   add!(mg, Demand(carrier = Electricity()),
                   Solar(),
                   Liion_rainflow(update_by_year = 12, soc_model = "linear", calendar = true,  fatigue_data = fatigue_data, couplage = (E=true, R=false)),
                   Grid(carrier = Electricity()))

   ω_d = Scenarios(mg, data["ω_optim"], adjust_length=true)

   i = 0
   for size_id in 1:Nb_sizing
      i+=1
      println(string(i,"/",Nb_sizing))
      for (couple_id, couple) in enumerate(couplages)
         println(couple)
         for (soc_id, SoC) in enumerate(SoC_model)
            println(SoC)
            for (soh_id, SoH) in enumerate(SoH_model)
               println(SoH)
               println("")


               #Select model and create microgrid
               mg = Microgrid(parameters = GlobalParameters(nh, ny, ns, renewable_share = 1.))

               if SoH == SoH_model[1]
                  #energie echangée
                  add!(mg, Demand(carrier = Electricity()),
                                  Solar(),
                                  Liion_energy_exchanged(;nCycle = 2 * fatigue_data.cycle[findfirst(fatigue_data.DoD .> 0.6)] , soc_model = SoC, couplage = couple, calendar = true),
                                  Grid(carrier = Electricity()))

               elseif SoH == SoH_model[2]
                  #rainflow

                  add!(mg, Demand(carrier = Electricity()),
                                  Solar(),
                                  Liion_rainflow(update_by_year = 12, soc_model = SoC, fatigue_data = fatigue_data, couplage = couple, calendar = true),
                                  Grid(carrier = Electricity()))

               elseif SoH == SoH_model[3]
                  #electrochimique
                  add!(mg, Demand(carrier = Electricity()),
                                  Solar(),
                                  Liion_electro_chimique(update_by_year = 12, soc_model = SoC, couplage = couple),
                                  Grid(carrier = Electricity()))

               else
                  add!(mg, Demand(carrier = Electricity()),
                                  Solar(),
                                  Liion_fixed_lifetime(soc_model = SoC, lifetime = 15),
                                  Grid(carrier = Electricity()))
               end

               #Simulation des microgrid créés

               #Le coût de l'abonnement va diluer les effets mais c'est la même pour tout le monde
               designer = initialize_designer!(mg, Manual(generations = [sizings[size_id,1]], storages = [sizings[size_id,2]], subscribed_power = [20]), ω_d)
               controller = initialize_controller!(mg, RBC(options = RBCOptions(policy_selection = 2)), ω_d)

               simulate!(mg, controller, designer, ω_d, options = Genesys.Options(mode = "serial", firstyear = false))

               metrics = Metrics(mg, designer)

               #Create ns*nb_metrics new rows and add them
               for scen in 1:ns

                     push!(value, sum(metrics.cost.total[:,scen]))
                     push!(metriques_names, "Cost")

                     push!(value, mean(metrics.renewable_share[2:ny,scen]))
                     push!(metriques_names, "RES")

                     push!(value,  sum(metrics.cost.total[:,scen]) - metrics.npv.salvage[ny,scen])
                     push!(metriques_names, "NPV*")

                     push!(value, metrics.npv.salvage[ny,scen])
                     push!(metriques_names, "Discounted salvage")

                     push!(value, sum(metrics.cost.opex[2:ny,scen]))
                     push!(metriques_names, "Discounted opex")

                     push!(value, sum(metrics.cost.capex[:,scen]))
                     push!(metriques_names, "Discounted capex")

                     push!(value, sum(metrics.npv.total[scen]))
                     push!(metriques_names, "NPV")

                     for tmp in 1:length(metriques_names_list)
                        push!(pv, sizings[size_id,1])
                        push!(bat, sizings[size_id,2])
                        push!(Couplage, couplages_names[couple_id])
                        push!(SOH, SoH_model_display[soh_id])
                        push!(SOC, SoC_model_display[soc_id])
                     end

               end


            end
         end
      end
   end


   df = Pandas.DataFrame(Dict("value" => value, "pv" => pv, "bat" => bat , "Couplage" => Couplage,  "SOH" => SOH, "SOC" => SOC, "Metrics" => metriques_names) )



   df_ref = df_pandas.pyo.get((df_pandas.SOH == SoH_model_display[3]) .& (df_pandas.SOC == SoC_model_display[2]) .& (df_pandas.Couplage == "ER"))
   value_ref = df_ref.value.to_numpy()




   for i in 1:length(SoH_model)
      for j in 1:length(SoC_model)
         figure(string(SoH_model[i] , " ", SoC_model[j]))
         df_tmp = Pandas.DataFrame(Dict("value" => [], "pv" => [], "bat" => [], "Couplage" => [],  "SOH" => [], "SOC" => [], "Metrics" => []) )

         for k in 1:length(couplages_names)



            df_filter = df.pyo.get((df.SOH == SoH_model_display[i]) .& (df.SOC == SoC_model_display[j]) .& (df.Couplage == couplages_names[k]))


            df_filter.insert(loc=1,
                column="Error (%)",
                value=(((df_filter.value - value_ref)/ value_ref).to_numpy()*100)
            )

            # df_filter.insert(loc=1,
            #    column="Flat Error",
            #    value=(((df_filter.value - value_ref)).to_numpy())
            # )

            #df_filter = df_filter.get((df_filter.Metrics != "NPV*") )



            if isempty(df_tmp)
               df_tmp = df_tmp.pyo[:append](df_filter, ignore_index = true )
            else
               df_tmp = df_tmp[:append](df_filter, ignore_index = true )
            end

         end
         ax = Seaborn.violinplot(x="Couplage", y="Error (%)", hue="Metrics", data=df_tmp, palette="Set2", inner="quartile")
         ax.axhline(0)
      end
   end



   df_pandas = Pandas.DataFrame(Dict("value" => value, "pv" => pv, "bat" => bat , "Couplage" => Couplage,  "SOH" => SOH, "SOC" => SOC, "Metrics" => metriques_names) )


   df_pandas = df_pandas.replace( "Linéaire", "Constant Yield")
   df_pandas = df_pandas.replace("Énergie échangée", "Exchanged Energy")
   df_pandas = df_pandas.replace( "Éléctro-chimique", "Electro-Chemical")
   df_pandas = df_pandas.replace( "Fixed", "Fixed Lifetime")

    SoC_model_display_EN = ["Constant Yield", "Tremblay-Dessaint"]
    SoH_model_display_EN = ["Exchanged Energy", "Rainflow", "Electro-Chemical", "Fixed Lifetime"]

   df_ref = df_pandas.get((df_pandas.SOH .== SoH_model_display_EN[3]) .& (df_pandas.SOC .== SoC_model_display_EN[2]) .& (df_pandas.Couplage .== "ER"))
   value_ref = df_ref.value.to_numpy()

#Keep only ER
   df_filter_ER = df_pandas.get(df_pandas.Couplage == couplages_names[4])
   df_filter_R = df_pandas.get(df_pandas.Couplage == couplages_names[3])
   df_filter_E = df_pandas.get(df_pandas.Couplage == couplages_names[2])
   df_filter = df_pandas.get(df_pandas.Couplage == couplages_names[1])



#Crée la metrique d'erreur

df_tmp = Pandas.DataFrame(Dict("value" => [], "pv" => [], "bat" => [], "Couplage" => [],  "SOC" => [], "SOH" => [], "Metrics" => []) )

   for i in 1:length(SoH_model)
      for j in 1:length(SoC_model)

            df_filter_soc_soh = df_filter_ER.get((df_filter_ER.SOH == SoH_model_display_EN[i]) .& (df_filter_ER.SOC == SoC_model_display_EN[j]))


            df_filter_soc_soh.insert(loc=1,
                column="Error (%)",
                value=(((df_filter_soc_soh.value - value_ref)/ value_ref).to_numpy()*100)
            )

            # df_filter.insert(loc=1,
            #    column="Flat Error",
            #    value=(((df_filter.value - value_ref)).to_numpy())
            # )

            #df_filter = df_filter.get((df_filter.Metrics != "NPV*") )



            if isempty(df_tmp)
               df_tmp = df_tmp.pyo[:append](df_filter_soc_soh, ignore_index = true )
            else
               df_tmp = df_tmp[:append](df_filter_soc_soh, ignore_index = true )
            end

      end
   end


#Garder uniquement les metriques NPV* et RES
   df_filter_ER = df_tmp.get((df_tmp.Metrics == "NPV*") .| (df_tmp.Metrics == "RES"))



   figure("Erreur en % des indicateurs techno-économiques par rapport à une référence")
   Seaborn.set_context("paper", rc=Dict("axes.labelsize"=>20, "legend.fontsize"=> 16,  "legend.title_fontsize"=> 20,
   "ytick.labelsize" => 16, "xtick.labelsize" => 16))

   ax  = Seaborn.catplot(x="SOH", y="Error (%)",
                   hue="Metrics", row="SOC",
                   data=df_filter_ER, kind="violin", split=true, inner="quartile",
                   palette="Set2",  order=["Fixed Lifetime","Exchanged Energy", "Rainflow", "Electro-Chemical"], bw=.2,
                   height=4, aspect=2.5);
   ax.set_xlabels("Ageing model")
   ax.map(plt.axhline, y=0, ls="-", c="black")


   # 2 différent figures

   df_filter_ER_lin = df_tmp.get(((df_tmp.Metrics == "NPV*") .| (df_tmp.Metrics == "RES")) .&  (df_tmp.SOC == "Constant Yield"))
   figure("Erreur en % des indicateurs techno-économiques par rapport à une référence")
   Seaborn.set_context("paper", rc=Dict("axes.labelsize"=>28, "legend.fontsize"=> 22,  "legend.title_fontsize"=> 28,
   "ytick.labelsize" => 22, "xtick.labelsize" => 22))

   ax  = Seaborn.catplot(x="SOH", y="Error (%)",
                   hue="Metrics",
                   data=df_filter_ER_lin, kind="violin", split=true, inner="quartile",
                   palette="Set2",  order=["Fixed Lifetime","Exchanged Energy", "Rainflow", "Electro-Chemical"], bw=.2,
                   height=4, aspect=2.5);
   ax.set_xlabels("Ageing model")
   ax.map(plt.axhline, y=0, ls="-", c="black")


   df_filter_ER_td = df_tmp.get(((df_tmp.Metrics == "NPV*") .| (df_tmp.Metrics == "RES")) .&  (df_tmp.SOC == "Tremblay-Dessaint"))
   figure("Erreur en % des indicateurs techno-économiques par rapport à une référence")
   Seaborn.set_context("paper", rc=Dict("axes.labelsize"=>28, "legend.fontsize"=> 22,  "legend.title_fontsize"=> 28,
   "ytick.labelsize" => 22, "xtick.labelsize" => 22))

   ax  = Seaborn.catplot(x="SOH", y="Error (%)",
                   hue="Metrics",
                   data=df_filter_ER_td, kind="violin", split=true, inner="quartile",
                   palette="Set2",  order=["Fixed Lifetime","Exchanged Energy", "Rainflow", "Electro-Chemical"], bw=.2,
                   height=4, aspect=2.5);
   ax.set_xlabels("Ageing model")
   ax.map(plt.axhline, y=0, ls="-", c="black")

#Le NPV total n'est pas représenté car son écart en % est potientiellement impertinent.
# cette métrique met en oeuvre des éléments de signes opposé et peut donc conduire à des résultat négatif, positif ou quasi nul,
# L'erreur ne devrait pas dans ce cas se mesurer en % d'erreur à la ref.

#Sauvegarde et loading des données.
using DataTables, TypedTables, IterableTables

new_df = DataFrames.DataFrame(tt)
CSV.write("Erreur_model_df2.csv", new_df)


df3 = CSV.read("Erreur_model_df2.csv", DataFrames.DataFrame)


replace!(df.SOC, "Linéaire" => "Constant Yield")
replace!(df.SOH, "Énergie échangée" => "Exchanged Energy")
replace!(df.SOH, "Éléctro-chimique" => "Electro-Chemical")
replace!(df.SOH, "Fixed" => "Fixed Lifetime")

tt = DataTable(df3)
df_pandas = Pandas.DataFrame(tt)
tt = DataTable(df_pandas)

df = DataFrames.DataFrame(tt)
