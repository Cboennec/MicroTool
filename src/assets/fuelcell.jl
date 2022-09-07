#=
    Fuel cell modelling
 =#

mutable struct FuelCell <: AbstractConverter
     # Paramètres
     α_p::Float64
     η_H2_E::Float64
     η_H2_H::Float64
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
     FuelCell(; α_p = 8/100,
             η_H2_E = 0.4,
             η_H2_H = 0.4,
             lifetime = 14,
             nHoursMax = 10000.,
             bounds = (lb = 0., ub = 50.),
             powerMax_ini = 0.,
             soh_ini = 1.) =
             new(α_p, η_H2_E, η_H2_H, lifetime, nHoursMax, bounds, powerMax_ini, soh_ini)
end

### Preallocation
function preallocate!(fc::FuelCell, nh::Int64, ny::Int64, ns::Int64)
     fc.powerMax = convert(SharedArray,zeros(ny+1, ns)) ; fc.powerMax[1,:] .= fc.powerMax_ini
     fc.carrier = [Electricity(), Heat(), Hydrogen()]
     fc.carrier[1].power = convert(SharedArray,zeros(nh, ny, ns))
     fc.carrier[2].power = convert(SharedArray,zeros(nh, ny, ns))
     fc.carrier[3].power = convert(SharedArray,zeros(nh, ny, ns))
     fc.soh = convert(SharedArray,zeros(nh+1, ny+1, ns)) ; fc.soh[1,1,:] .= fc.soh_ini
     fc.cost = convert(SharedArray,zeros(ny, ns))
     return fc
end

### Operation dynamic
function compute_operation_dynamics!(h::Int64, y::Int64, s::Int64, fc::FuelCell, decision::Float64, Δh::Int64)
    fc.soh[h+1,y,s], fc.carrier[1].power[h,y,s], fc.carrier[2].power[h,y,s], fc.carrier[3].power[h,y,s] =
    compute_operation_dynamics(fc, (powerMax = fc.powerMax[y,s], soh = fc.soh[h,y,s]), decision, Δh)
end

function compute_operation_dynamics(fc::FuelCell, state::NamedTuple{(:powerMax, :soh), Tuple{Float64, Float64}}, decision::Float64, Δh::Int64)
    # Power constraint and correction
    fc.α_p * state.powerMax <= decision && state.soh * fc.nHoursMax / Δh > 1. ? power_E = min(decision, state.powerMax) : power_E = 0.
    # Power conversion
    power_H2 = - power_E / fc.η_H2_E
    power_H = - power_H2 * fc.η_H2_H
    # SoH computation
    soh_next = state.soh - (power_E > 0.) * Δh / fc.nHoursMax
    return soh_next, power_E, power_H, power_H2
end

### Investment dynamic
function compute_investment_dynamics!(y::Int64, s::Int64, fc::FuelCell, decision::Union{Float64, Int64})
    fc.powerMax[y+1,s], fc.soh[1,y+1,s] = compute_investment_dynamics(fc, (powerMax = fc.powerMax[y,s], soh = fc.soh[end,y,s]), decision)
end

function compute_investment_dynamics(fc::FuelCell, state::NamedTuple{(:powerMax, :soh), Tuple{Float64, Float64}}, decision::Union{Float64, Int64})
    if decision > 1e-2
        powerMax_next = decision
        soh_next = 1.
    else
        powerMax_next = state.powerMax
        soh_next = state.soh
    end
    return powerMax_next, soh_next
end
