#=
    Thermal energy storage modelling
 =#

mutable struct ThermalStorage <: AbstractStorage
     # Paramètres
     α_p_ch::Float64
     α_p_dch::Float64
     η_ch::Float64
     η_dch::Float64
     η_self::Float64
     α_soc_min::Float64
     α_soc_max::Float64
     lifetime::Int64
     bounds::NamedTuple{(:lb, :ub), Tuple{Float64, Float64}}
     # Initial conditions
     Erated_ini::Float64
     soc_ini::Float64
     soh_ini::Float64
     # Variables
     Erated::AbstractArray{Float64,2}
     carrier::Heat
     soc::AbstractArray{Float64,3}
     # Eco
     cost::AbstractArray{Float64,2}
     # Inner constructor
     ThermalStorage(; α_p_ch = 1.5,
          α_p_dch = 1.5,
          η_ch = 0.8,
          η_dch = 0.8,
          η_self = 0.008,
          α_soc_min = 0.,
          α_soc_max = 1.,
          lifetime = 25,
          bounds = (lb = 0., ub = 1000.),
          Erated_ini = 1e-6,
          soc_ini = 0.5,
          soh_ini = 1.) =
          new(α_p_ch, α_p_dch, η_ch, η_dch, η_self, α_soc_min, α_soc_max, lifetime, bounds, Erated_ini, soc_ini, soh_ini)
end

### Preallocation
function preallocate!(tes::ThermalStorage, nh::Int64, ny::Int64, ns::Int64)
    tes.Erated = convert(SharedArray,zeros(ny+1, ns)) ; tes.Erated[1,:] .= tes.Erated_ini
    tes.carrier = Heat()
    tes.carrier.power = convert(SharedArray,zeros(nh, ny, ns))
    tes.soc = convert(SharedArray,zeros(nh+1, ny+1, ns)) ; tes.soc[1,1,:] .= tes.soc_ini
    tes.cost = convert(SharedArray,zeros(ny, ns))
    return tes
end

### Operation dynamic
function compute_operation_dynamics!(h::Int64, y::Int64, s::Int64, tes::ThermalStorage, decision::Float64, Δh::Int64)
     tes.soc[h+1,y,s], tes.carrier.power[h,y,s] = compute_operation_dynamics(tes, (Erated = tes.Erated[y,s], soc = tes.soc[h,y,s]), decision, Δh)
end

function compute_operation_dynamics(tes::ThermalStorage, state::NamedTuple{(:Erated, :soc), Tuple{Float64, Float64}}, decision::Float64, Δh::Int64)
     # Control power constraint and correction
     power_dch = max(min(decision, tes.α_p_dch * state.Erated, tes.η_dch * (state.soc * (1. - tes.η_self * Δh) - tes.α_soc_min) * state.Erated / Δh), 0.)
     power_ch = min(max(decision, -tes.α_p_ch * state.Erated, (state.soc * (1. - tes.η_self * Δh) - tes.α_soc_max) * state.Erated / Δh / tes.η_ch), 0.)
     # SoC dynamic
     soc_next = state.soc * (1. - tes.η_self * Δh) - (power_ch * tes.η_ch + power_dch / tes.η_dch) * Δh / state.Erated
     return soc_next, power_dch + power_ch
end

### Investment dynamic
function compute_investment_dynamics!(y::Int64, s::Int64, tes::ThermalStorage, decision::Union{Float64, Int64})
    tes.Erated[y+1,s], tes.soc[1,y+1,s] = compute_investment_dynamics(tes, (Erated = tes.Erated[y,s], soc = tes.soc[end,y,s]), decision)
end

function compute_investment_dynamics(tes::ThermalStorage, state::NamedTuple{(:Erated, :soc), Tuple{Float64, Float64}}, decision::Union{Float64, Int64})
    if decision > 1e-2
        Erated_next = decision
        soc_next = tes.soc_ini
    else
        Erated_next = state.Erated
        soc_next = state.soc
    end
    return Erated_next, soc_next
end
