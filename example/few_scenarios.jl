
# Load packages
include("..\\..\\src\\Genesys.jl")

using Main.Genesys
using JLD, Dates, Seaborn
using CSV, DataFrames
using Statistics

pygui(true)

# Parameters of the simulation
const nh, ny, ns = 8760, 16, 5

# Load input data
data = load(joinpath("example","data","ausgrid_5_twostage.jld"))

data = load(joinpath("data_light.jld"))

using HDF5, JLD


new_pv = (t = data2["pv"].t[:,:,1:10], power = data2["pv"].power[:,:,1:10], cost = data2["pv"].cost[:,1:10])
new_elyz = (cost =  data2["elyz"].cost[:,1:10],)
new_liion = (cost = data2["liion"].cost[:,1:10],)
new_fc = (cost =  data2["fc"].cost[:,1:10],)
new_tes = (cost = data2["tes"].cost[:,1:10],)
new_heater = (cost =  data2["heater"].cost[:,1:10],)
new_h2tank = (cost = data2["h2tank"].cost[:,1:10],)
new_grid = (cost_in = data2["grid"].cost_in[:,:,1:10], cost_out = data2["grid"].cost_out[:,:,1:10])
new_ld_H = (t = data2["ld_H"].t[:,:,1:10], power = data2["ld_H"].power[:,:,1:10])
new_ld_E = (t = data2["ld_E"].t[:,:,1:10], power = data2["ld_E"].power[:,:,1:10])

new_data = Dict("pv"=>new_pv, "elyz"=>new_elyz, "liion"=>new_liion, "fc"=>new_fc, "tes"=>new_tes, "heater"=>new_heater, "h2tank"=>new_h2tank, "grid"=>new_grid, "ld_H"=>new_ld_H, "ld_E"=>new_ld_E)

save("data_light.jld", "data", new_data)
