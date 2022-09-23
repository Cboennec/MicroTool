
# Load packages
include("..\\..\\src\\Genesys.jl")

using Main.Genesys
using JLD, Dates, Seaborn
using CSV, DataFrames
using Profile
using Combinatorics

pygui(true)

# Parameters of the simulation
const nh, ny, ns = 8760, 20, 100

# Load input data
data = load(joinpath("example","data","ausgrid_5_twostage.jld"))

#https://www.researchgate.net/figure/Battery-cycle-to-failure-versus-the-DOD-gap-34_fig5_276136726
fatigue_data = DataFrame(CSV.File("example\\data\\fatigue_data.csv.csv", delim = ";", header = [Symbol("DoD"),Symbol("cycle")], types=Dict(:DoD=>Float64, :cycle=>Float64)))



# Initialize the microgrid
microgrid_vermeer = Microgrid(parameters = GlobalParameters(nh, ny, ns, renewable_share = 1.))
microgrid_linear = Microgrid(parameters = GlobalParameters(nh, ny, ns, renewable_share = 1.))
microgrid_rainflow = Microgrid(parameters = GlobalParameters(nh, ny, ns, renewable_share = 1.))
microgrid_electrochimique = Microgrid(parameters = GlobalParameters(nh, ny, ns, renewable_share = 1.))

#Add the equipment to the microgrid
add!(microgrid_vermeer, Demand(carrier = Electricity()),
                Solar(),
                Liion_vermeer(soc_model = "vermeer"),
                Grid(carrier = Electricity()))

add!(microgrid_linear, Demand(carrier = Electricity()),
                Solar(),
                Liion_energy_exchanged(soc_model = "linear"),
                Grid(carrier = Electricity()))

add!(microgrid_rainflow, Demand(carrier = Electricity()),
                Solar(),
                Liion_rainflow(soc_model = "tremblay_dessaint"),
                Grid(carrier = Electricity()))

add!(microgrid_electrochimique, Demand(carrier = Electricity()),
                Solar(),
                Liion_electro_chimique(soc_model = "tremblay_dessaint"),
                Grid(carrier = Electricity()))

# Initialize scenarios
 ω_a = Scenarios(microgrid_vermeer, data["ω_simu"], adjust_length=true)
 ω_b = Scenarios(microgrid_linear, data["ω_simu"], adjust_length=true)
 ω_c = Scenarios(microgrid_rainflow, data["ω_simu"], adjust_length=true)
 ω_d = Scenarios(microgrid_electrochimique, data["ω_simu"], adjust_length=true)


# Initialize the designer
controller = RBC(options = RBCOptions(policy_selection =  2))
designer = Manual(generations = [20.], storages = [20.], subscribed_power = [5.])

designer_vermeer = initialize_designer!(microgrid_vermeer, designer, ω_a)
controller_vermeer = initialize_controller!(microgrid_vermeer, controller, ω_a)
#
designer_linear = initialize_designer!(microgrid_linear, designer, ω_b)
controller_linear = initialize_controller!(microgrid_linear, controller, ω_b)

designer_rainflow = initialize_designer!(microgrid_rainflow, designer, ω_c)
controller_rainflow = initialize_controller!(microgrid_rainflow, controller, ω_c)

designer_electrochimique = initialize_designer!(microgrid_electrochimique, designer, ω_d)
controller_electrochimique = initialize_controller!(microgrid_electrochimique, controller, ω_d)

set = [1,2,3,4]
perms = collect(permutations(set))

time1, time2, time3, time4 = 0,0,0,0

@time simulate!(microgrid_rainflow, controller_rainflow, designer_rainflow, ω_c, options = Genesys.Options(mode = "multithreads"))
@time simulate!(microgrid_electrochimique, controller_electrochimique, designer_electrochimique, ω_d, options = Genesys.Options(mode = "multithreads"))
@time simulate!(microgrid_vermeer, controller_vermeer, designer_vermeer, ω_a, options = Genesys.Options(mode = "multithreads"))
@time simulate!(microgrid_linear, controller_linear, designer_linear, ω_b, options = Genesys.Options(mode = "multithreads"))


for perm in perms
    for id in perm
        if id == 1
          time_ver += (@timed simulate!(microgrid_vermeer, controller_vermeer, designer_vermeer, ω_a, options = Genesys.Options(mode = "multithreads")) ).time
        end
        if id == 2
          time_lin += (@timed simulate!(microgrid_linear, controller_linear, designer_linear, ω_b, options = Genesys.Options(mode = "multithreads")) ).time
        end
        if id == 3
          time_rain += (@timed simulate!(microgrid_rainflow, controller_rainflow, designer_rainflow, ω_c, options = Genesys.Options(mode = "multithreads")) ).time
        end
        if id == 4
          time_elec += (@timed simulate!(microgrid_electrochimique, controller_electrochimique, designer_electrochimique, ω_d, options = Genesys.Options(mode = "multithreads")) ).time
        end

    end
end

nb_perms = length(perms)

time_rain /= nb_perms
time_elec /= nb_perms
time_ver /= nb_perms
time_lin /= nb_perms


#
# open("prof_rain.txt", "w") do s
#     Profile.print(IOContext(s, :displaysize => (24, 500)))
# end
# Juno.profiler()

 # Compute the metrics
#
# open("prof_ele.txt", "w") do s
#     Profile.print(IOContext(s, :displaysize => (24, 500)))
# end
# Juno.profiler()


#
# open("prof_ver.txt", "w") do s
#     Profile.print(IOContext(s, :displaysize => (24, 500)))
# end
# Juno.profiler()


#
# open("prof_lin.txt", "w") do s
#     Profile.print(IOContext(s, :displaysize => (24, 500)))
# end
# Juno.profiler()


metrics = Metrics(microgrid_electrochimique, designer_electrochimique)

# Plots
plot_operation(microgrid_electrochimique, y = 2:ny, s=1)
plot_metrics(metrics)
