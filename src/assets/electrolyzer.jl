#=
    Electrolyzer modelling
 =#

mutable struct Electrolyzer <: AbstractConverter
     # Paramètres
     α_p::Float64
     η_E_H2::Float64
     η_E_H::Float64
     lifetime::Int64
     nHoursMax::Float64
     bounds::NamedTuple{(:lb, :ub), Tuple{Float64, Float64}}
     # Initial conditions
     powerMax_ini::Float64
     soh_ini::Float64
     # Variables
     powerMax::AbstractArray{Float64,2}
     carrier::Vector{EnergyCarrier}
     soh::AbstractArray{Float64,3}
     # Eco
     cost::AbstractArray{Float64,2}
     # Inner constructor
     Electrolyzer(; α_p = 5/100,
                 η_E_H2 = 0.5,
                 η_E_H = 0.3,
                 lifetime = 15,
                 nHoursMax = 26000.,
                 bounds = (lb = 0., ub = 50.),
                 powerMax_ini = 0.,
                 soh_ini = 1.) =
                 new(α_p, η_E_H2, η_E_H, lifetime, nHoursMax, bounds, powerMax_ini, soh_ini)
end

### Preallocation
function preallocate!(elyz::Electrolyzer, nh::Int64, ny::Int64, ns::Int64)
     elyz.powerMax = convert(SharedArray,zeros(ny+1, ns)) ; elyz.powerMax[1,:] .= elyz.powerMax_ini
     elyz.carrier = [Electricity(), Heat(), Hydrogen()]
     elyz.carrier[1].power = convert(SharedArray,zeros(nh, ny, ns))
     elyz.carrier[2].power = convert(SharedArray,zeros(nh, ny, ns))
     elyz.carrier[3].power = convert(SharedArray,zeros(nh, ny, ns))
     elyz.soh = convert(SharedArray,zeros(nh+1, ny+1, ns)) ; elyz.soh[1,1,:] .= elyz.soh_ini
     elyz.cost = convert(SharedArray,zeros(ny, ns))
     return elyz
end

### Operation dynamic
function compute_operation_dynamics!(h::Int64, y::Int64, s::Int64, elyz::Electrolyzer, decision::Float64, Δh::Int64)
    elyz.soh[h+1,y,s], elyz.carrier[1].power[h,y,s], elyz.carrier[2].power[h,y,s], elyz.carrier[3].power[h,y,s] =
    compute_operation_dynamics(elyz, (powerMax = elyz.powerMax[y,s], soh = elyz.soh[h,y,s]), decision, Δh)
end

function compute_operation_dynamics(elyz::Electrolyzer, state::NamedTuple{(:powerMax, :soh), Tuple{Float64, Float64}}, decision::Float64, Δh::Int64)
    # Power constraint and correction
    elyz.α_p * state.powerMax >= decision && state.soh * elyz.nHoursMax / Δh > 1. ? power_E = max(decision, -state.powerMax) : power_E = 0.
    # Power conversion
    power_H2 = - power_E * elyz.η_E_H2
    power_H = - power_E * elyz.η_E_H
    # SoH computation
    soh_next = state.soh - (power_E > 0.) * Δh / elyz.nHoursMax
    return soh_next, power_E, power_H, power_H2
end
### Investment dynamic
function compute_investment_dynamics!(y::Int64, s::Int64, elyz::Electrolyzer, decision::Union{Float64, Int64})
    elyz.powerMax[y+1,s], elyz.soh[1,y+1,s] = compute_investment_dynamics(elyz, (powerMax = elyz.powerMax[y,s], soh = elyz.soh[end,y,s]), decision)
end

function compute_investment_dynamics(elyz::Electrolyzer, state::NamedTuple{(:powerMax, :soh), Tuple{Float64, Float64}}, decision::Union{Float64, Int64})
    if decision > 1e-2
        powerMax_next = decision
        soh_next = 1.
    else
        powerMax_next = state.powerMax
        soh_next = state.soh
    end
    return powerMax_next, soh_next
end
