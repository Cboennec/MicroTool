
# Load packages
include("..\\..\\src\\Genesys.jl")

using Main.Genesys
using JLD, Dates, Seaborn
using CSV, DataFrames, Pandas

pygui(true)

# Parameters of the simulation
const nh, ny, ns = 8760, 20, 100

# Load input data
data = load(joinpath("example","data","ausgrid_5_twostage.jld"))

#https://www.researchgate.net/figure/Battery-cycle-to-failure-versus-the-DOD-gap-34_fig5_276136726
fatigue_data = DataFrame(CSV.File("example\\data\\fatigue_data.csv.csv", delim = ";", header = [Symbol("DoD"),Symbol("cycle")], types=Dict(:DoD=>Float64, :cycle=>Float64)))



# Initialize the microgrid
microgrid_linear = Microgrid(parameters = GlobalParameters(nh, ny, ns, renewable_share = 1.))
microgrid_vermeer = Microgrid(parameters = GlobalParameters(nh, ny, ns, renewable_share = 1.))


#Add the equipment to the microgrid

add!(microgrid_linear, Demand(carrier = Electricity()),
                Solar(),
                Liion_energy_exchanged(soc_model = "linear"),
                Grid(carrier = Electricity()))
add!(microgrid_vermeer, Demand(carrier = Electricity()),
                Solar(),
                Liion_vermeer(soc_model = "vermeer"),
                Grid(carrier = Electricity()))
# Initialize scenarios
ω_a = Scenarios(microgrid_vermeer, data["ω_simu"], adjust_length=true)
ω_b = Scenarios(microgrid_linear, data["ω_simu"], adjust_length=true)


nh = 8760

hours = []
scenarios = []
PVs = []
loads = []
for h in 2:nh
        hours = vcat(hours,repeat([h],ns))
        scenarios = vcat(scenarios, [a for a in 1:ns])
        loads = vcat(loads,ω_a.demands[1].power[h,2,:])
        PVs = vcat(PVs,ω_a.generations[1].power[h,2,:])
end

df = Pandas.DataFrame(Dict(:Hours=>hours, :scenario=>scenarios, :Load=>loads, :PV=>PVs))


plotLoad = Seaborn.lineplot(data=df, x="Hours", y="Load")
plt.xlabel("Hours")
plt.ylabel("Load [kWh]")
plotLoad.get_children()[1].set_color("red")


plotPV = Seaborn.lineplot(data=df, x="Hours", y="PV")
plt.xlabel("Hours")
plt.ylabel("PV")
plotLoad.get_children()[1].set_color("red")






nh = 8760

hours = []
scenarios = []
type = []
value = []
for data_type in 1:2
        for h in 1:nh
                hours = vcat(hours,repeat([h],ns))
                scenarios = vcat(scenarios, [a for a in 1:ns])
                if data_type == 1
                        value = vcat(value,ω_a.demands[1].power[h,2,:])
                        type = vcat(type,repeat(["Load"],ns))
                else
                        value = vcat(value,ω_a.generations[1].power[h,2,:])
                        type = vcat(type,repeat(["PV"],ns))
                end
        end
end

df = Pandas.DataFrame(Dict(:Hours=>hours, :scenario=>scenarios, :type=>type, :value=>value))


g = Seaborn.FacetGrid(df, row="type",
                  height=1.7, aspect=4,)
g.map(Seaborn.lineplot,"Hours", "value")
Seaborn.despine()
