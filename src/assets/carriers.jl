abstract type EnergyCarrier end

mutable struct Electricity <: EnergyCarrier
    power::AbstractArray{Float64,3}
    Electricity() = new()
end

mutable struct Heat <: EnergyCarrier
    power::AbstractArray{Float64,3}
    Heat() = new()
end

mutable struct Hydrogen <: EnergyCarrier
    power::AbstractArray{Float64,3}
    Hydrogen() = new()
end
