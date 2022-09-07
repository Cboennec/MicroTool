module Genesys
# TODO: break into submodules
# Optimisation
using JuMP, Cbc, Metaheuristics, SDDP
# Math
using Statistics, StatsBase, MultivariateStats, Clustering, Distributions, Distances, LinearAlgebra, Interpolations
# Others
using Seaborn, ProgressMeter, Dates, Distributed, SharedArrays, CSV, DataFrames, JLD
# Assets
include(joinpath("assets","microgrid.jl"))
include(joinpath("assets","carriers.jl"))
include(joinpath("assets","liion","liion.jl"))
include(joinpath("assets","liion","liion_energy_exchanged.jl"))
include(joinpath("assets","liion","liion_rainflow.jl"))
include(joinpath("assets","liion","liion_electro_chimique.jl"))
include(joinpath("assets","liion","liion_vermeer.jl"))
include(joinpath("assets","liion","liion_fixed_lifetime.jl"))
include(joinpath("assets","tes.jl"))
include(joinpath("assets","h2tank.jl"))
include(joinpath("assets","electrolyzer.jl"))
include(joinpath("assets","fuelcell.jl"))
include(joinpath("assets","heater.jl"))
include(joinpath("assets","grid.jl"))
include(joinpath("assets","solar.jl"))
include(joinpath("assets","demand.jl"))
export AbstractController, AbstractLiion, AbstractDesigner
export Microgrid, Demand, Solar, Liion_energy_exchanged, Liion_rainflow, Liion_fixed_lifetime, Liion_vermeer, Liion_electro_chimique, Tremblay_dessaint_params, vermeer_params, Electro_chimique_params, ThermalStorage, H2Tank, FuelCell, Electrolyzer, Heater, Grid, GlobalParameters
export Electricity, Heat, Hydrogen
export add!
# Scenarios
include(joinpath("scenarios","scenarios.jl"))
include(joinpath("scenarios","reduction.jl"))
include(joinpath("scenarios","generation.jl"))
include(joinpath("scenarios","utils.jl"))
export Scenarios
export ManualReducer, SAAReducer, MeanValueReducer, FeatureBasedReducer
export UnitRangeTransform, ZScoreTransform
export PCAReduction, StatsReduction
export KmedoidsClustering
export MarkovGenerator, AnticipativeGenerator
export reduce, generate
# Optimization utils
include(joinpath("optimization","utils.jl"))
export Expectation, CVaR, WorstCase
# Operation optimization
include(joinpath("optimization","controller","dummy.jl"))
include(joinpath("optimization","controller","anticipative.jl"))
include(joinpath("optimization","controller","rb.jl"))
include(joinpath("optimization","controller","olfc.jl"))
export Dummy, RBC, Anticipative, OLFC
export RBCOptions, AnticipativeOptions, OLFCOptions
export initialize_controller!
# Investment optimization
include(joinpath("optimization","designer","manual.jl"))
include(joinpath("optimization","designer","milp.jl"))
include(joinpath("optimization","designer","metaheuristic.jl"))
export Manual, Metaheuristic, MILP
export MetaheuristicOptions, MILPOptions, ManualOptions
export initialize_designer!
# Simulation
include(joinpath("simulation","informations.jl"))
include(joinpath("simulation","dynamics.jl"))
include(joinpath("simulation","power_balances.jl"))
include(joinpath("simulation","simulations.jl"))
export simulate!
# Utils
include(joinpath("utils","metrics.jl"))
include(joinpath("utils","plots.jl"))
include(joinpath("utils","saves.jl"))
export Metrics
export plot_operation, plot_metrics
export COST

end
